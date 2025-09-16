#!/usr/bin/env python

from collections import defaultdict
import csv
import sys


class KrakenEntry:
    """
    A class representing a line in a kraken report.

    Attributes:
        taxon_id (str): The NCBI taxon identifier.
        name (string): The scientific name associated with this taxon.
        rank (str): A letter coding the rank of this taxon.
        depth (int): The number of indentations this entry had in the kraken file (related to hierarchy).
        count (int): The count of reads assigned to this taxon and its descendants.
        ucount (int): The count of reads assigned specifically to this taxon.
        domain (str): The domain this taxon is a member of (name not taxon_id).
        parent (str): The taxon id associated with the taxonomic parent.
        children (set): A set of taxon ids associated with the direct taxonomic children.
        sibling_rank (int): An integer representing the ranking among direct siblings (share the parent) based on count.
        hierarchy (list): An ordered list of taxon ids representing the parents to taxonomic root.
    """

    def __init__(self, row=None, domain=None, hierarchy=[]):
        """
        Initializes an KrakenEntry object.

        Parameters:
            row (str): A row from a kraken file
            domain (str): The taxonomic domain this entry is associated with.
            hierarchy (list): A list of taxon_ids tracing the ancestry from this taxon to the root.
        """
        self.taxon_id = "0"
        self.name = "unclassified"
        self.rank = "U"
        self.depth = 0
        self.count = 0  # inclusive count
        self.ucount = 0  # unique_count
        self.domain = domain
        self.parent = None
        self.children = set()
        self.sibling_rank = 0
        self.hierarchy = hierarchy
        if row is not None:
            self.add_row(row)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def print(self):
        """
        Print the attributes of KrakenEntry as a string
        """
        print(
            f"{self.taxon_id},{self.name},{self.rank},{self.depth},{self.count},{self.ucount},{self.domain},{self.parent},{self.children},{self.sibling_rank},{self.hierarchy}"
        )

    def parse_depth(self, name):
        """
        Parse the number of indentations from the kraken name string

        Args:
            name (str): The "Scientific Name" column from a kraken report line

        Returns:
            depth (int): The relative depth within the name string.
        """
        parse_name = name.split(" ")
        depth = 0
        for i in parse_name:
            if i != "":
                break
            depth += 1
        depth = int(depth / 2)
        return depth

    def add_row(self, row):
        """
        Parse information from kraken report row to update this entry.

        Args:
            row (str): A line from a kraken report.
        """
        self.taxon_id = row["Taxonomy ID"]
        self.name = row["Scientific Name"].strip()
        self.depth = self.parse_depth(row["Scientific Name"])
        self.rank = row["Rank"]
        self.count = int(row["Clades"])  # inclusive count
        self.ucount = int(row["Taxonomies"])  # unique_count
        if self.count < self.ucount:
            self.count, self.ucount = self.ucount, self.count
        self.hierarchy = self.hierarchy[: self.depth]

    def add_parent(self, parent):
        assert self.parent == None or parent == self.parent
        self.parent = parent

    def add_child(self, child):
        self.children.add(child)

    def set_sibling_rank(self, rank):
        self.sibling_rank = rank

    def update(self, new_entry):
        """
        Checks if the information in an existing (this) entry matches the new entry.
        Does not change the counts or ucounts.
        Allows name, rank and children to be updated.

        Args:
            new_entry (KrakenEntry): A KrakenEntry object.
        """
        if self.name != new_entry.name:
            print(f"Updated name {self.name} to {new_entry.name}")
            self.name = new_entry.name
        if self.rank != new_entry.rank:
            print(f"Updated rank {self.rank} to {new_entry.rank}")
            self.rank = new_entry.rank
        assert self.depth == new_entry.depth

        # self.count += new_entry.count
        # self.ucount += new_entry.ucount

        assert self.domain == new_entry.domain
        assert self.parent == new_entry.parent
        self.children.update(new_entry.children)

        assert self.hierarchy == new_entry.hierarchy


class KrakenReport:
    """
    A class representing a kraken report.

    Attributes:
        entries (dict): A dict with keys for taxon ids and values for KrakenEntry.
        total (int): Total number of reads in report.
        unclassified (int): Number of unclassified reads.
        classified (int): Number of classified reads.
        domains (int): A dict with keys for names of domains and values for associated taxon id.
        file_name (Path): File path for report
    """

    def __init__(self, file_name=None):
        """
        Initializes an KrakenReport object.

        Parameters:
            file_name (Path): Name of kraken report file.
        """
        self.entries = defaultdict(KrakenEntry)
        self.total = 0
        self.unclassified = 0
        self.classified = 0
        self.domains = defaultdict(str)  # maps name to taxon_id
        self.file_name = file_name
        if file_name:
            self.load_file(file_name)
            self.unclassified = self.entries["0"].count
            self.classified = self.entries["1"].count if "1" in self.entries else 0
            self.total = self.classified + self.unclassified

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def print(self):
        """
        Print the attributes of KrakenEntry as a string
        """
        print(
            f"Report has {len(self.entries)} taxon entries corresponding to {self.classified} classified and {self.unclassified} unclassified reads."
        )

    def add_parent_child(self, parent_id, child_id):
        """
        Add parent-child relationship to each KrakenEntry.

        Parameters:
            parent_id (str): taxon_id of parent.
            child_id (str): taxon_id of child.
        """

        self.entries[child_id].add_parent(parent_id)
        self.entries[parent_id].add_child(child_id)

    def set_sibling_ranks(self):
        """
        Rank siblings (share common parent) based on the number of classified reads (count including descendants). Lower
        rank means higher read count. Rank starts at 1, 2, 3, ...
        Only operates below domain level

        Parameters:
            parent_id (str): taxon_id of parent.
            child_id (str): taxon_id of child.
        """
        for entry_id, entry in self.entries.items():
            if entry.sibling_rank > 0 or entry.parent is None:
                continue
            if entry.rank in ["D", "R", "R1"]:
                entry.set_sibling_rank(1)
            elif len(self.entries[entry.parent].children) == 1:
                entry.set_sibling_rank(1)
            else:
                sibling_dict = {
                    i: self.entries[i].count
                    for i in self.entries[entry.parent].children
                }
                sorted_counts = sorted(sibling_dict.values(), reverse=True)
                for i, c in sibling_dict.items():
                    rank = sorted_counts.index(c) + 1
                    self.entries[i].set_sibling_rank(rank)

    def check_sibling_ranks(self):
        """
        Check every entry has been set a sibling_rank. 0 means not set.
        """
        for entry_id in self.entries:
            if self.entries[entry_id].sibling_rank == 0:
                print(entry_id)
                assert entry_id in ["0", "1"]

    def check_report(self, file_name):
        """
        Check every line in the kraken report has the appropriate number of tab separated fields

        Parameters:
            file_name (Path): Name of kraken report file

        Returns:
            report_has_header (bool): True if report includes the header line
            num_fields (int) number of fields in a line [6,8]
        """
        with open(file_name, "r") as handle:
            line = handle.readline()
            num_fields = len(line.split("\t"))
            if num_fields not in [6, 8]:
                sys.stderr.write(
                    f"Kraken report file {file_name} badly formatted - must have 6 or 8 columns"
                )
                sys.exit(9)
            if line.startswith("%"):
                return True, num_fields
            else:
                return False, num_fields

    def load_file(self, file_name):
        """
        Check every line in the kraken report has the appropriate number of tab separated fields

        Parameters:
            file_name (Path): Name of kraken report file

        Returns:
            report_has_header (bool): True if report includes the header line
            num_fields (int) number of fields in a line [6,8]
        """
        csvfile = open(file_name, newline="")
        df = {}
        report_has_header, num_fields = self.check_report(file_name)
        if report_has_header:
            df = csv.DictReader(csvfile, delimiter="\t")
        elif num_fields == 6:
            df = csv.DictReader(
                csvfile,
                delimiter="\t",
                fieldnames=[
                    "% of Seqs",
                    "Clades",
                    "Taxonomies",
                    "Rank",
                    "Taxonomy ID",
                    "Scientific Name",
                ],
            )
        elif num_fields == 8:
            df = csv.DictReader(
                csvfile,
                delimiter="\t",
                fieldnames=[
                    "% of Seqs",
                    "Clades",
                    "Taxonomies",
                    "Read Minimizers",
                    "Taxon Minimizers",
                    "Rank",
                    "Taxonomy ID",
                    "Scientific Name",
                ],
            )

        hierarchy = []
        domain = None
        for row in df:
            try:
                if row["Rank"] == "D":
                    domain = row["Scientific Name"].strip()
                    self.domains[domain] = row["Taxonomy ID"]
                entry = KrakenEntry(row=row, domain=domain, hierarchy=hierarchy)

            except:
                sys.stderr.write(
                    f"Found badly formatted row:\n{row}\n. Quitting load of {file_name}."
                )
                sys.exit(9)

            self.entries[entry.taxon_id] = entry
            hierarchy = entry.hierarchy.copy()
            if len(hierarchy) > 0:
                self.add_parent_child(hierarchy[-1], entry.taxon_id)
            if entry.taxon_id != "0":
                hierarchy.append(entry.taxon_id)
        csvfile.close()
        self.set_sibling_ranks()
        # self.check_sibling_ranks()

    def get_domains(self):
        """
        Get a list of taxonomic domains found in the report

        Returns:
            list: List of domains
        """
        domains = []
        for entry_id, entry in self.entries.items():
            if entry.rank == "D":
                domains.append(entry_id)
                entry.print()
        return domains

    def get_tips(self):
        """
        Get a list of terminating taxa (ie have no children)

        Returns:
            list: List of KrakenEntry
        """
        tips = []
        for entry_id, entry in self.entries.items():
            if len(entry.children) == 0 and entry_id != "0":
                tips.append(entry_id)
                entry.print()
        return tips

    def get_rank_entries(self, rank):
        """
        Get a list of  taxa at a given rank

        Parameters:
            rank (str): a taxonomic rank
        Returns:
            list: List of KrakenEntry
        """
        subset = []
        for entry_id, entry in self.entries.items():
            if entry.rank == rank:
                subset.append(entry_id)
                entry.print()
        return subset

    def get_percentage(self, taxon_id, denominator="classified"):
        """
        For a taxon_id, what percentage of the dataset counts (or the domain counts) is it?

        Parameters:
            taxon_id (str): a taxon ID
            denominator (str): Optional name of a domain or ["classified", "total"]
        Returns:
            float: Percentage of the dataset corresponding to counts of this taxon and descendants
        """
        if denominator == "classified":
            total = self.classified
        elif denominator == "total":
            total = self.total
        elif denominator in self.domains:
            total = self.entries[self.domains[denominator]].count
        else:
            print(f"Not a valid denominator {denominator}")

        if (
            denominator not in ["classified", "total"]
            and self.entries[taxon_id].domain != denominator
        ):
            return 0.0

        count = self.entries[taxon_id].count
        percentage = float(count) / float(total) * 100
        # print(f"{count}/{total} = {percentage:.2f}")
        return percentage

    def to_source_target_df(
        self, out_file="source_target.csv", max_rank=None, domain=None
    ):
        """
        Convert this KrakenReport to a CSV with "source", "target", "value", "percentage" columns. If max_rank given,
        the result is filtered to including only this number of top-ranking children for each parent. If domain is
        given, includes only those as part of this domain, and scales percentage accordingly. Output file written to
        "source_target.csv".

        Parameters:
            max_rank (str): (optional) maximum rank among siblings to be included
            domain (str): (optional) name of a domain
        Returns:
            list: a list of row-entry dictionaries
        """
        records = []
        ignore = set()
        skip = set()
        denominator = "total"

        for entry_id, entry in self.entries.items():

            # we don't want the connections to cellular organisms, root etc
            if entry.parent in [None, "0", "1", "131567"]:
                continue

            # filter by domain where required
            if domain:
                denominator = domain
                if entry.domain != domain:
                    continue

            # filter by rank when specified
            if max_rank:
                if entry.sibling_rank > max_rank:
                    ignore.add(entry_id)
                    continue
                elif entry.parent in ignore:
                    ignore.add(entry_id)
                    continue

            # filter if an intermediate rank
            if entry.rank not in [
                "K",
                "D",
                "D1",
                "D2",
                "P",
                "C",
                "O",
                "F",
                "G",
                "G1",
                "S",
                "S1",
                "S2",
            ]:
                skip.add(entry_id)
                continue

            index = 1
            while index < len(entry.hierarchy) and entry.hierarchy[-index] in skip:
                index += 1
            source_id = entry.hierarchy[-index]
            records.append(
                {
                    "source": self.entries[source_id].name,
                    "target": entry.name,
                    "value": entry.count,
                    "percentage": self.get_percentage(
                        entry_id, denominator=denominator
                    ),
                }
            )

        with open(f"{out_file}", "w", newline="") as csvfile:
            fieldnames = ["source", "target", "value", "percentage"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in records:
                writer.writerow(row)

        print(len(records), len(ignore), len(skip))

        return records

    def to_df(self, sample_id="sample_id", ranks=[]):
        """
        Create a pandas dataframe from KrakenReport. If ranks are specified, include only entries at those ranks.

        Parameters:
            sample_id (str): (optional) A name to give this report in the dataframe.
            ranks (list): (optional) A list of strings corresponding to ranks to include.
        """
        import pandas as pd

        if not ranks or len(ranks) == 0:
            taxon_ids = [e for e in self.entries.keys()]
        else:
            taxon_ids = [
                e for e in self.entries.keys() if self.entries[e].rank in ranks
            ]
        return pd.DataFrame(
            {sample_id: [self.entries[e].count for e in taxon_ids]}, index=taxon_ids
        )

    def check_host(self, host_dict):
        """
        Check if any host taxon_id exceeds the max number of reads permitted.

        Parameters:
            host_dict (dict): A dictionary with key a taxon_id, value a max number of reads allowed.
        """
        for host_id, max_host_count in host_dict.items():
            if host_id in self.entries and self.entries[host_id].count > max_host_count:
                sys.stderr.write(
                    f"ERROR: found {self.entries[host_id].count} reads corresponding to host {self.entries[host_id].name} with taxon_id {host_id}, max allowed is {max_host_count}\n"
                )
                sys.exit(2)

    def update_entry(self, new_entry):
        """
        Adds a KrakenEntry to the entries dictionary, checking if it already exists that the contents match.
        Note that the result DOES NOT UPDATE/TRANSFER COUNTS

        Parameters:
            new_entry (KrakenEntry): A new kraken entry object.
        """
        # print(f"updating {new_entry.taxon_id}")
        if new_entry.taxon_id in self.entries:
            self.entries[new_entry.taxon_id].update(new_entry)
        else:
            self.entries[new_entry.taxon_id] = new_entry
            self.entries[new_entry.taxon_id].ucount = 0
            self.entries[new_entry.taxon_id].count = 0

    def get_mrca(self, taxon_id_1, taxon_id_2):
        """
        Find the mrca taxon_id from the hierarchy lists between 2 taxon_ids in the entries dictionary.

        Parameters:
            taxon_id_1 (str): First taxon_id in entries.
            taxon_id_2 (str): Second taxon_id in entries.

        Returns:
            mrca_taxon_id
        """
        i = 0
        assert taxon_id_1 in self.entries and taxon_id_2 in self.entries
        if taxon_id_1 == "0":
            # print(f"MRCA of old {taxon_id_1} and new {taxon_id_2} is {taxon_id_1}")
            return taxon_id_1

        entry1 = self.entries[taxon_id_1]
        entry2 = self.entries[taxon_id_2]

        if taxon_id_1 == "1" and taxon_id_1 in entry2.hierarchy:
        #    print(f"MRCA of old {taxon_id_1} and new {taxon_id_2} is {taxon_id_1}")
            return taxon_id_1

        while (
            i < len(entry1.hierarchy)
            and i < len(entry2.hierarchy)
            and entry1.hierarchy[i] == entry2.hierarchy[i]
        ):
            if i == len(entry1.hierarchy) - 1 or i == len(entry2.hierarchy) - 1:
                break
            elif entry1.hierarchy[i+1] != entry2.hierarchy[i+1]:
                break
            i += 1

        # print(f"MRCA of old {taxon_id_1} and new {taxon_id_2} is {entry1.hierarchy} position {i}")
        return entry1.hierarchy[i]

    def update_counts(self, changes):
        """
        Uses a dictionary of changes to update the counts and ucounts

        Parameters:
            changes (dict): A dictionary mapping old_taxon_id, new_taxon_id to number of counts transferred from old to new.
        """
        for old_taxon_id in changes:
            for new_taxon_id in changes[old_taxon_id]:
                print(f"Moving {changes[old_taxon_id][new_taxon_id]} counts from {old_taxon_id} to {new_taxon_id}")
                
                mrca = self.get_mrca(old_taxon_id, new_taxon_id)
                print(f"MRCA of {old_taxon_id} and {new_taxon_id} is {mrca}")

                self.entries[old_taxon_id].ucount -= changes[old_taxon_id][new_taxon_id]
                print(f"Removing {changes[old_taxon_id][new_taxon_id]} ucounts from {old_taxon_id}")

                if not (old_taxon_id == "1" and old_taxon_id in self.entries[new_taxon_id].hierarchy):
                    self.entries[old_taxon_id].count -= changes[old_taxon_id][new_taxon_id]
                    print(f"Removing {changes[old_taxon_id][new_taxon_id]} counts from {old_taxon_id}")

                assert self.entries[old_taxon_id].ucount >= 0

                if old_taxon_id != "0":
                    for taxon_id in reversed(self.entries[old_taxon_id].hierarchy):
                        if taxon_id != mrca:
                            self.entries[taxon_id].count -= changes[old_taxon_id][
                                new_taxon_id
                            ]
                            print(f"Removing {changes[old_taxon_id][new_taxon_id]} counts from {taxon_id}")
                            assert self.entries[taxon_id].count >= 0
                        elif taxon_id == mrca:
                            break

                if not (new_taxon_id == "1" and new_taxon_id in self.entries[old_taxon_id].hierarchy):
                    self.entries[new_taxon_id].ucount += changes[old_taxon_id][new_taxon_id]
                    self.entries[new_taxon_id].count += changes[old_taxon_id][new_taxon_id]
                    print(f"Adding {changes[old_taxon_id][new_taxon_id]} counts and ucounts to {new_taxon_id}")

                for taxon_id in reversed(self.entries[new_taxon_id].hierarchy):
                    if taxon_id != mrca:
                        self.entries[taxon_id].count += changes[old_taxon_id][
                            new_taxon_id
                        ]
                        print(f"Adding {changes[old_taxon_id][new_taxon_id]} counts to {taxon_id}")

                    elif taxon_id == mrca:
                        break

                self.unclassified = self.entries["0"].count
                self.classified = self.entries["1"].count if "1" in self.entries else 0
                if (self.total != self.classified + self.unclassified):
                    print(f"Broke after {old_taxon_id} and {new_taxon_id} with {self.unclassified}, {self.classified}")
                    assert self.total == self.classified + self.unclassified

    def clean(self):
        """
        Removes entries which have 0 counts and references to them.
        """
        set_zeroes = set()
        for taxon_id in self.entries:
            if self.entries[taxon_id].count == 0:
                set_zeroes.add(taxon_id)
        for taxon_id in set_zeroes:
            entry = self.entries[taxon_id]
            if entry.parent in self.entries:
                print(taxon_id, entry.parent, self.entries[entry.parent].children)
                self.entries[entry.parent].children.remove(taxon_id)
            for child in entry.children:
                if child in self.entries:
                    assert child in set_zeroes
            del self.entries[taxon_id]

    def update(self, new_report, changes):
        """
        This process updates the existing KrakenReport with entries from the new KrakenReport.
        If the existing one is empty, entries are just copied over. Otherwise...
        new (zero count) KrakenEntry objects are added into the existing entries
        structure (and common KrakenEntry objects are checked for compatibility with name, rank
        and children updated where required). Counts are changed based on a prescribed dictionary
        of `changes`. Any zero count entries are removed and sibling ranks reevaluated.
        Total, unclassified and classified are updated.

        Parameters:
            new_report (KrakenReport): A new loaded kraken report.
        """
        print(f"New report has {new_report.entries.keys()} keys")
        if len(self.entries) == 0:
            self.entries = new_report.entries
            self.unclassified = self.entries["0"].count
            self.classified = self.entries["1"].count if "1" in self.entries else 0
            self.total = self.classified + self.unclassified
            print(f"Merged report has {self.entries.keys()} keys")
        else:
            print(
                f"New report has {len(new_report.entries)} items and existing report has {len(self.entries)} items"
            )
            for taxon_id, new_entry in new_report.entries.items():
                self.update_entry(new_entry)
            self.update_counts(changes)
            self.clean()
            self.set_sibling_ranks()
            self.unclassified = self.entries["0"].count
            self.classified = self.entries["1"].count if "1" in self.entries else 0
            assert self.total == self.classified + self.unclassified

    def save(self, file_name=None):
        """
        Save the KrakenReport object in kraken report format
        """
        if not file_name:
            file_name = self.file_name
        with open(file_name, "w") as out:
            fieldnames = [
                "% of Seqs",
                "Clades",
                "Taxonomies",
                "Rank",
                "Taxonomy ID",
                "Scientific Name",
            ]
            header = "\t".join(fieldnames)
            out.write(f"{header}\n")
            for taxon_id, entry in self.entries.items():
                name_buff = "  " * entry.depth
                percentage = self.get_percentage(taxon_id, "total")
                line = f"{percentage:6.2f}\t{entry.count}\t{entry.ucount}\t{entry.rank}\t{entry.taxon_id}\t{name_buff}{entry.name}"
                out.write(f"{line}\n")
