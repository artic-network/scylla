#!/usr/bin/env python
"""
Single source of truth for the taxonomic rank vocabulary used across the report.

NCBI restructured its top-level ranks in 2024: `superkingdom` was retired, with
Bacteria/Archaea/Eukaryota becoming rank `domain` and Viruses becoming
`acellular root` (with `realm` clades such as Duplodnaviria beneath it). Both the
old and new names are listed here so a lineage tops out at a real domain under
either an old or a current taxonomy dump.

This module is stdlib-free so it can be imported by scripts running in minimal
task containers. It is consumed by:
  - aggregate_lineages_bracken.py, to decide which lineage nodes to keep;
  - report.py, for the kraken-report domain-name heuristic;
  - make_report.py, which injects RANKS and DOMAIN_RANKS into the report's
    JavaScript so the two sides cannot drift apart.
"""

# Ordered coarse -> fine. `cellular root` (taxid 131567, "cellular organisms") is
# deliberately absent: including it would collapse Bacteria, Archaea and
# Eukaryota under a single meaningless top-level node. Rank `no rank` (taxid 1,
# "root") is likewise absent.
RANKS = [
    "acellular root",
    "superkingdom",
    "domain",
    "realm",
    "clade",
    "kingdom",
    "phylum",
    "subphylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "subspecies",
    "serotype",
]

# The ranks that mark a domain boundary, i.e. the top level of a lineage. A node
# at one of these ranks names the domain its descendants belong to.
DOMAIN_RANKS = ["acellular root", "superkingdom", "domain"]

# The rank reported for the synthetic "Unclassified" node. It is rendered as a
# top-level node alongside the real domains, so it must be one of DOMAIN_RANKS
# for the report's domain grouping to pick it up.
UNCLASSIFIED_RANK = "domain"

# Domain/realm names used to identify the domain boundary in a *kraken report*,
# which carries single-letter rank codes rather than the rank strings above. This
# is a name-based heuristic and is intentionally separate from DOMAIN_RANKS: see
# the comment at its use site in report.py.
KNOWN_DOMAIN_NAMES = {"Bacteria", "Archaea", "Eukaryota", "Viruses"}
