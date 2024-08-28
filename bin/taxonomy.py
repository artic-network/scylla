import os
import sys
from collections import defaultdict


class Taxonomy:
    """
    A class representing taxonomic information.

    Attributes:
        parents (dict): A dict with keys for taxon ids and values for parent taxon id.
        children (dict): A dict with keys for taxon ids and values for sets of direct child taxon id.
    """
    def __init__(self, taxonomy_dir):
        self.parents = {}
        self.children = defaultdict(set)

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
                    tax_id, parent_tax_id = fields[0], fields[1]

                    self.parents[tax_id] = parent_tax_id
                    self.children[parent_tax_id].add(tax_id)
        except:
            sys.stderr.write(
                "ERROR: Could not find taxonomy nodes.dmp file in %s" % taxonomy_dir
            )
            sys.exit(4)
