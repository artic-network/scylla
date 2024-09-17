import sys
from collections import defaultdict


def trim_read_id(read_id):
    if read_id.endswith("/1") or read_id.endswith("/2"):
        read_id = read_id[:-2]

    return read_id


def parse_kraken_assignment_line(line):
    line_vals = line.strip().split("\t")
    if len(line_vals) < 5:
        return -1, ""
    if "taxid" in line_vals[2]:
        temp = line_vals[2].split("taxid ")[-1]
        taxon_id = temp[:-1]
    else:
        taxon_id = line_vals[2]

    read_id = trim_read_id(line_vals[1])

    if taxon_id == "A":
        taxon_id = "81077"
    else:
        taxon_id = taxon_id
    return taxon_id, read_id


class KrakenAssignments:
    """
    A class representing a kraken assignment file.

    Attributes:
        file (str): Name of file to parse.
    """
    def __init__(self, assignment_file):
        self.file = assignment_file

    def parse_kraken_assignment_file(self, taxon_id_map, parents=None):
        """
        Parses the kraken assignment file and collects the read_ids associated with each of the
        required taxon ids.

        Parameters:
            taxon_id_map (dict): A dict from a taxon id to a list of related taxon_ids.
            parents (dict): A dict mapping taxon id to parent taxon id from NCBI Taxonomy

        Returns:
            read_map (dict): A dict from read_id to a set of values from the taxon_id_map if the
                             read_id was classified as the taxon_id.
        """
        read_map = defaultdict(set)
        with open(self.file, "r") as kfile:
            for line in kfile:
                taxon_id, read_id = parse_kraken_assignment_line(line)
                if taxon_id in taxon_id_map:
                    read_map[read_id] = {taxon_id} #.update(taxon_id_map[taxon_id])
                elif parents:
                    # handle case where taxon_id has changed
                    current = taxon_id
                    while current in parents and current not in taxon_id_map and current != "1":
                        current = parents[current]
                        if current in taxon_id_map:
                            print(f"Add {taxon_id} to {current} list")
                            read_map[read_id] = current #.update(taxon_id_map[current])
        return read_map

