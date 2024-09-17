import os
import sys
from collections import defaultdict

class TaxonEntry:
    """
    A class representing a line in a kraken report.

    Attributes:
        taxon_id (str): The NCBI taxon identifier.
        name (string): The scientific name associated with this taxon.
        rank (str): A letter coding the rank of this taxon.
    """
    def __init__(self, taxon_id="0", name="unclassified", rank="U"):
        self.taxon_id = taxon_id
        self.name = name
        self.rank = rank

    def print(self):
        print(
            f"{self.taxon_id},{self.name},{self.rank}")


class Taxonomy:
    """
    A class representing taxonomic information.

    Attributes:
        parents (dict): A dict with keys for taxon ids and values for parent taxon id.
        children (dict): A dict with keys for taxon ids and values for sets of direct child taxon id.
    """
    def __init__(self, taxonomy_dir=None, taxon_ids=None):
        self.parents = defaultdict(str)
        self.children = defaultdict(set)
        self.entries = defaultdict(TaxonEntry)

        if taxonomy_dir:
            self.load_parents_and_children(taxonomy_dir)
        if taxonomy_dir and taxon_ids:
            self.load_entries(taxonomy_dir, taxon_ids)

    def load_parents_and_children(self, taxonomy_dir):
        """
        Loads the parent child relationships from the "nodes.dmp" file in the taxonomy directory.
        Updates the class members `parents` and `children`

        Parameters:
            taxonomy_dir (str): The unzipped directory downloaded from NCBI taxonomy.
        """
        taxonomy = os.path.join(taxonomy_dir, "nodes.dmp")
        try:
            with open(taxonomy, "r") as f:
                for line in f:
                    fields = line.split("\t|\t")
                    taxon_id, parent_taxon_id = fields[0], fields[1]
                    self.parents[taxon_id] = parent_taxon_id
                    self.children[parent_taxon_id].add(taxon_id)
        except:
            sys.stderr.write(
                f"ERROR: Could not find taxonomy nodes.dmp file in {taxonomy_dir}"
            )
            sys.exit(4)

    def load_entries_from_nodes(self, taxonomy_dir, taxon_ids):
        if len(taxon_ids) == 0:
            return

        taxonomy = os.path.join(taxonomy_dir, "nodes.dmp")
        if not os.path.exists(taxonomy):
            sys.stderr.write(
                f"ERROR: Could not find taxonomy nodes.dmp file in {taxonomy_dir}"
            )
            sys.exit(4)
        try:
            with open(taxonomy, "r") as f:
                for line in f:
                    fields = line.split("\t|\t")
                    taxon_id, rank = fields[0], fields[2]
                    self.entries[taxon_id].taxon_id = taxon_id
                    self.entries[taxon_id].rank = rank
        except:
            sys.stderr.write(
                f"ERROR: Badly formatted nodes.dmp file in {taxonomy}"
            )
            sys.exit(4)

    def load_entries_from_names(self, taxonomy_dir, taxon_ids):
        if len(taxon_ids) == 0:
            return

        taxonomy = os.path.join(taxonomy_dir, "names.dmp")
        if not os.path.exists(taxonomy):
            sys.stderr.write(
                f"ERROR: Could not find taxonomy names.dmp file in {taxonomy_dir}"
            )
            sys.exit(4)
        try:
            with open(taxonomy, "r") as f:
                for line in f:
                    fields = [i.lstrip() for i in line.split("\t|")]
                    taxon_id, name, name_type = fields[0], fields[1], fields[3]
                    if taxon_id in taxon_ids and ("scientific name" in name_type or self.entries[taxon_id].name == "unclassified"):
                        self.entries[taxon_id].name = name
        except:
            sys.stderr.write(
                f"ERROR: Badly formatted names.dmp file in {taxonomy}"
            )
            sys.exit(4)

    def load_entries(self, taxonomy_dir, taxon_ids):
        self.load_entries_from_nodes(taxonomy_dir, taxon_ids)
        self.load_entries_from_names(taxonomy_dir, taxon_ids)

    def get_taxon_id_map(self, taxon_ids=[], include_unclassified=False):
        """
        Generates a map from a specified list of taxon_ids and their children to sets of taxon_ids
        Any descendent of a specified taxon_id will include this taxon_id in it's set.
        If include_unclassified is specified, there will be an entry from 0 to the input taxon_ids

        Parameters:
            taxonomy_dir (str): The unzipped directory downloaded from NCBI taxonomy.
            taxon_ids (list/set): An iterable of taxon_ids to consider.
        """
        taxon_id_map = defaultdict(set)

        for key in taxon_ids:
            taxon_id_map[key].add(key)
        if include_unclassified:
            taxon_id_map["0"].update(taxon_ids)

        check = list(taxon_ids)
        while len(check) > 0:
            current = check.pop()
            check.extend(self.children[current])
            for child in self.children[current]:
                taxon_id_map[child].update(taxon_id_map[current])

        return taxon_id_map