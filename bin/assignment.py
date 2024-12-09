#!/usr/bin/env python

import sys
from collections import defaultdict


def trim_read_id(read_id):
    """
    Parses the read_id to remove forward/reverse identifier.

    Parameters:
        read_id (str): A read name.

    Returns:
        read_id (str): Trimmed read name without forward/reverse identifiers.
    """
    if read_id.endswith("/1") or read_id.endswith("/2"):
        read_id = read_id[:-2]

    return read_id

def get_mrca(taxon_id1, taxon_id2, parents):
    if taxon_id1 == taxon_id2:
        return taxon_id1
    elif taxon_id1 == "0" or taxon_id2 == "0":
        return "0"

    ancestry1 = []
    current = taxon_id1
    while current in parents and current != "1":
        ancestry1.append(current)
        current = parents[current]
    ancestry1.append(current)
    ancestry1.reverse()

    ancestry2 = []
    current = taxon_id2
    while current in parents and current != "1":
        ancestry2.append(current)
        current = parents[current]
    ancestry2.append(current)
    ancestry2.reverse()

    index = min(len(ancestry1), len(ancestry2)) - 1
    while index > 0 and ancestry1[index] != ancestry2[index]:
        index -= 1
    return ancestry1[index]


class KrakenAssignmentEntry:
    """
    A class representing a line in a kraken assignment file.

    Attributes:
        classified (str): C if read classified, U if unclassified
        read_id (str): The read name
        taxon_id (str): The NCBI taxon identifier
        length (int): Length of read in bp
        kmer_string (str): space separated string representing the taxon_ids matched along the read
    """
    def __init__(self, line=None):
        """
        Initializes an KrakenAssignmentEntry object.

        Parameters:
            row (str): A row from a kraken file
        """
        self.classified = "U"
        self.read_id = ""
        self.taxon_id = "0"
        self.length = 0
        self.kmer_string = ""
        if line is not None:
            self.add_line(line)
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def add_line(self, line):
        """
        Parses the line in the kraken assignment file.

        Parameters:
            line (str): A line from kraken assignment file.

        """
        num_fields = len(line.split("\t"))
        if num_fields != 5:
            sys.stderr.write(
                f"Kraken assignment line {line} badly formatted - must have 5 fields"
            )
            sys.exit(11)
        self.classified, self.read_id, self.taxon_id, length, self.kmer_string = line.strip().split("\t")
        self.length = int(length)

        #if "taxid" in self.taxon_id:
        #    temp = self.taxon_id.split("taxid ")[-1]
        #    self.taxon_id = temp[:-1] // can't remember where this came from so leave it out

        self.read_id = trim_read_id(self.read_id)

        if self.taxon_id == "A":
            self.taxon_id = "81077"

    def declassify(self):
        """
            Changes the classified status of KrakenAssignmentEntry to unclassified
        """
        self.classified = "U"
        self.taxon_id = "0"

    def get_line(self):
        """
            Get string representation of KrakenAssignmentEntry
        """
        fields = [self.classified, self.read_id, self.taxon_id, str(self.length), self.kmer_string]
        return "\t".join(fields)

    def print(self):
        """
        Print the attributes of KrakenAssignmentEntry as a string
        """
        print(f"{self.get_line()}")


class KrakenAssignments:
    """
    A class representing a kraken assignment file.

    Attributes:
        file_name (str): Name of file to parse.
        load (bool): If set loads the contents of the file into memory
    """
    def __init__(self, assignment_file, load=False):
        """
        Initializes an KrakenAssignments object.
        """
        self.entries = defaultdict(KrakenAssignmentEntry)
        self.file_name = assignment_file

        if load:
            self.load_file()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def get_read_map(self, taxon_id_map, parents=None):
        """
        Parses the kraken assignment file and collects the read_ids associated with each of the
        required taxon ids. If paired reads are provided, will consider the common ancestor of
        the 2 assignments

        Parameters:
            taxon_id_map (iter): Iterable of taxon ids to identify reads for.
            parents (dict): A dict mapping taxon id to parent taxon id from NCBI Taxonomy

        Returns:
            read_map (dict): A dict from read_id to a taxon_id in the input iterable.
        """
        read_map = defaultdict(str)
        extended_map = defaultdict(str)
        comments = set()
        with open(self.file_name, "r") as kfile:
            for line in kfile:
                assignment = KrakenAssignmentEntry(line)
                taxon_id, read_id = assignment.taxon_id, assignment.read_id

                corrected_taxon_id = taxon_id
                if parents:
                    while corrected_taxon_id in parents and corrected_taxon_id not in taxon_id_map and corrected_taxon_id != "1":
                        corrected_taxon_id = parents[corrected_taxon_id]
                        if corrected_taxon_id in taxon_id_map and corrected_taxon_id!= taxon_id:
                            comments.add(f"Assign {taxon_id} to {corrected_taxon_id} list")

                if read_id in extended_map:
                    mrca_taxon_id = get_mrca(corrected_taxon_id, extended_map[read_id], parents)
                    if mrca_taxon_id in taxon_id_map:
                        if mrca_taxon_id != read_map[read_id]:
                            comments.add(f"Reassign {extended_map[read_id]} (and {corrected_taxon_id}) to mrca {mrca_taxon_id} list")
                            read_map[read_id] = mrca_taxon_id
                    elif read_id in read_map:
                        comments.add(f"MRCA {mrca_taxon_id} of {extended_map[read_id]} and {corrected_taxon_id} not in taxon_id_map")
                        del read_map[read_id]

                elif corrected_taxon_id in taxon_id_map:
                    if corrected_taxon_id != taxon_id:
                        comments.add(f"Assign {taxon_id} to {corrected_taxon_id} list")
                    read_map[read_id] = corrected_taxon_id
                extended_map[read_id] = taxon_id
        for c in comments:
            print(c)
        return read_map

    def load_file(self, taxon_ids=None):
        """
        Loads all entries in the kraken assignment file. If this is a paired file and there is a clash, result is unclassified

        Parameters:
            taxon_ids (iterable): A subset of taxon_ids to retain assignment lines from.
        """
        with open(self.file_name, "r") as kfile:
            for line in kfile:
                assignment = KrakenAssignmentEntry(line)
                if (taxon_ids and assignment.taxon_id in taxon_ids) or not taxon_ids:
                    if assignment.read_id in self.entries and assignment.taxon_id != self.entries[assignment.read_id].taxon_id:
                        self.entries[assignment.read_id].declassify()
                    else:
                        self.entries[assignment.read_id] = assignment

    def update(self, new_assignments, changes=None):
        """
        Updates read assignments using new file with preference

        Parameters:
            new_assignments (KrakenAssignments): A new loaded KrakenAssignments object.
        """
        if not changes:
            changes = defaultdict(lambda : defaultdict(int))

        if len(self.entries) == 0:
            self.entries = new_assignments.entries
            return

        for read_id, entry in new_assignments.entries.items():
            if read_id not in self.entries:
                self.entries[read_id] = entry
                changes["0"][entry.taxon_id] += 1

            elif (read_id in self.entries
                  and entry.classified == "C"
                  and entry.taxon_id != self.entries[read_id].taxon_id):
                old_taxon_id = self.entries[read_id].taxon_id
                new_taxon_id = entry.taxon_id
                self.entries[read_id] = entry
                changes[old_taxon_id][new_taxon_id] += 1

        return changes

    def save(self):
        """
            Save the KrakenAssignments object in kraken assignment format
        """
        with open(self.file_name, "w") as out:
            for taxon_id, entry in self.entries.items():
                out.write(f"{entry.get_line()}\n")

