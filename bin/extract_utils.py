#!/usr/bin/env python

import sys
import gzip
import pyfastx
import json
from collections import defaultdict, OrderedDict
from pathlib import Path
import statistics as stats

from assignment import trim_read_id

import numpy as np


# Credit to @nrminor from the μbioinfo slack for this approach!
def ascii_to_phred(q: str, offset: int = 33) -> np.ndarray:
    # q: ASCII quality string from a single FASTQ record
    # offset: 33 for Sanger / modern Illumina, 64 for ancient Illumina
    arr = np.frombuffer(q.encode("ascii"), dtype=np.uint8)
    return arr - offset


def mean_phred(q: str, offset: int = 33) -> float:
    """
    Mean Phred quality for a single FASTQ quality string. Unlike np.mean(ascii_to_phred(q)),
    this never materialises a full subtracted array - summing the raw ASCII byte values
    (numpy upcasts the accumulator, so this doesn't overflow) and subtracting the offset once
    from the mean is the same value, computed in about a third of the time.
    """
    if not q:
        return 0.0
    arr = np.frombuffer(q.encode("ascii"), dtype=np.uint8)
    return float(arr.sum()) / arr.size - offset


def record_to_fastq_entry(record):
    """
    Converts a record to a FASTQ entry string.

    Args:
        record (tuple): A tuple containing (name, seq, qual).
    """
    name, seq, qual = record
    return f"@{name}\n{seq}\n+\n{qual}\n"


def check_read_files(reads):
    """
    Takes a read filename and checks if the file is zipped, a FASTA or FASTQ format.

    Args:
        reads (str): The reads filename.

    Returns:
        filetype (str): Either "fasta" or "fastq".
        zipped (bool): Whether the file is gzipped or not.
    """
    if reads[-3:] == ".gz":
        read_file = gzip.open(reads, "rt")
        zipped = True
    else:
        read_file = open(reads, "rt")
        zipped = False
    first = read_file.readline()
    if len(first) == 0:
        sys.stderr.write("ERROR: sequence file's first line is blank\n")
        sys.exit(5)
    if first[0] == ">":
        filetype = "fasta"
    elif first[0] == "@":
        filetype = "fastq"
    else:
        sys.stderr.write("ERROR: sequence file must be FASTA or FASTQ\n")
        sys.exit(5)
    return filetype, zipped


def setup_outfiles(fastq_2, prefixes, filetype):
    """
    Sets up the output filenames for each key (usually a taxon ID).

    Args:
        fastq_2 (bool): True or a string if fastq_2 exists, False or None if reads were not paired.
        prefixes (dict): A dictionary from the key (usually taxon_id) to the prefix of the output file for that key.
        filetype (str): "fasta" or "fastq" depending on input filetype.

    Returns:
        filenames (dict): Dictionary with keys from the prefixes, to a list of output filenames.
    """
    filenames = defaultdict(list)

    for key, prefix in prefixes.items():
        if fastq_2:
            filenames[key].append(f"{prefix}_1.{filetype}")
            filenames[key].append(f"{prefix}_2.{filetype}")
        else:
            filenames[key].append(f"{prefix}.{filetype}")
    return filenames


def _default_max_open_files():
    try:
        import resource

        soft, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
        return max(16, min(128, soft // 4))
    except Exception:
        return 64


class TaxonWriter:
    """
    Buffers extracted FASTQ/FASTA records per taxon ID in memory and periodically flushes them to disk. Output file
    handles are pooled with a bounded LRU cache so the number of simultaneously open files never exceeds `max_open`,
    regardless of how many distinct taxon IDs are being written to - this avoids hitting the OS open-file-descriptor
    limit while also avoiding buffering every read for every taxon in memory until the end of the file (as previously
    happened for the low-memory=False code path).

    Records for a given taxon are always flushed in the order they were written, so read order within each output
    file is preserved exactly as before.
    """

    def __init__(
        self,
        filenames,
        file_index,
        max_open=None,
        per_taxon_bytes=1 << 20,
        total_bytes=512 << 20,
    ):
        self.filenames = filenames
        self.file_index = file_index
        self.max_open = max_open if max_open else _default_max_open_files()
        self.per_taxon_bytes = per_taxon_bytes
        self.total_bytes = total_bytes

        self.buffers = defaultdict(list)
        self.buffer_bytes = defaultdict(int)
        self.total_buffered = 0
        self.handles = OrderedDict()  # taxon_id -> open file handle, LRU order

    def write(self, taxon_id, record):
        entry = record_to_fastq_entry(record)
        self.buffers[taxon_id].append(entry)
        n = len(entry)
        self.buffer_bytes[taxon_id] += n
        self.total_buffered += n

        if self.buffer_bytes[taxon_id] >= self.per_taxon_bytes:
            self._flush(taxon_id)
        elif self.total_buffered >= self.total_bytes:
            self._flush_largest()

    def _flush_largest(self):
        while self.total_buffered >= self.total_bytes and self.buffer_bytes:
            taxon_id = max(self.buffer_bytes, key=self.buffer_bytes.get)
            self._flush(taxon_id)

    def _flush(self, taxon_id):
        if not self.buffers[taxon_id]:
            return
        handle = self._handle(taxon_id)
        handle.write("".join(self.buffers[taxon_id]))
        self.total_buffered -= self.buffer_bytes[taxon_id]
        self.buffers[taxon_id] = []
        self.buffer_bytes[taxon_id] = 0

    def _handle(self, taxon_id):
        if taxon_id in self.handles:
            self.handles.move_to_end(taxon_id)
            return self.handles[taxon_id]

        if len(self.handles) >= self.max_open:
            _, evicted = self.handles.popitem(last=False)
            evicted.close()

        # Append mode is safe even for the first write to a given path: each
        # writer operates on a fresh Nextflow task working directory, so the
        # output file never pre-exists.
        path = self.filenames[taxon_id][self.file_index]
        handle = open(path, "a")
        self.handles[taxon_id] = handle
        return handle

    def close(self):
        for taxon_id in list(self.buffers.keys()):
            self._flush(taxon_id)
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()


def new_taxon_stats():
    """
    A fresh, empty running-stats accumulator for a single taxon.
    """
    return {"count": 0, "sum_qual": 0.0, "sum_len": 0}


def update_summary_with_record(taxon_id, record, out_counts, read_stats):
    """
    Updates the summary dictionaries for counts, qualities and lengths of reads from the new record. Rather than
    keeping every individual read's quality/length, a running count/sum per taxon is kept - this bounds memory use
    regardless of how many reads are seen for a given taxon.

    Args:
        taxon_id (str): Key for dictionaries
        record (SeqRecord): A read parsed by Bio.SeqIO
        out_counts (dict): Taxon ID to read count for that taxon.
        read_stats (dict): Taxon ID to a running {"count", "sum_qual", "sum_len"} accumulator for that taxon.
    """
    out_counts[taxon_id] += 1
    name, seq, qual = record

    taxon_stats = read_stats[taxon_id]
    taxon_stats["count"] += 1
    taxon_stats["sum_qual"] += mean_phred(qual)
    taxon_stats["sum_len"] += len(seq)


def file_iterator(
    read_file,
    read_map,
    subtaxa_map,
    inverse,
    out_counts,
    read_stats,
    writer,
):
    """
    Iterate through the read_file file and add the read to the appropriate output file, via the writer, updating
    summary statistics as it goes.

    Args:
        read_file (str): Name of read file.
        read_map (dict): Dictionary from read ID to list of taxon ID associated with it.
        subtaxa_map (dict): Dictionary from taxon ID (output by read map) to list of taxon ID associated with
                            out files to be updated
        inverse (bool): If True, exclude all reads which match the taxon IDs in the read_map. If False, select them.
        out_counts (dict): Taxon ID to read count for that taxon.
        read_stats (dict): Taxon ID to running {"count", "sum_qual", "sum_len"} accumulator for that taxon.
        writer (TaxonWriter): Buffered, bounded-filehandle writer for this file (forward or reverse read).
    Returns:
        int: Count of reads written to file.
    """
    count = 0
    total_length = 0

    sys.stderr.write(f"Reading in {read_file}\n")
    for record in pyfastx.Fastq(read_file, build_index=False):
        name, seq, qual = record
        total_length += len(seq)
        trimmed_name = trim_read_id(name)
        # A single dict.get() replaces the old `set(read_map.keys())` membership
        # check followed by a second read_map[...] lookup - read_map is a
        # defaultdict(str), so .get() (which never triggers the default
        # factory) is used rather than [] to distinguish "not present" (None)
        # from a stored value.
        taxon_id = read_map.get(trimmed_name)
        if inverse:
            if taxon_id is not None:
                for taxon in subtaxa_map[taxon_id]:
                    update_summary_with_record(taxon, record, out_counts, read_stats)
                continue
            else:
                count += 1
                update_summary_with_record("other", record, out_counts, read_stats)
                writer.write("other", record)

        elif not inverse:
            if taxon_id is None:
                continue
            for taxon in subtaxa_map[taxon_id]:
                count += 1
                update_summary_with_record(taxon, record, out_counts, read_stats)
                writer.write(taxon, record)

    writer.close()

    return count, total_length


def file_iterator_multi(read_file, contexts):
    """
    Like file_iterator, but applies several independent fraction contexts to a single pass over
    read_file, instead of requiring one pass per context. Each context is a dict with the same
    per-fraction state file_iterator would otherwise take as separate arguments: 'read_map',
    'subtaxa_map', 'inverse', 'out_counts', 'read_stats', 'writer'. Every context's inclusion/
    exclusion logic is applied to every record exactly as file_iterator would for that context
    alone - contexts never interact (each has its own read_map/writer/stats), only the expensive
    per-record decompress+parse work is shared.

    Args:
        read_file (str): Name of read file.
        contexts (list): List of per-fraction context dicts, as above.
    Returns:
        counts (list): Count of reads written to file, one per context, in the same order.
        total_length (int): Total length of all reads in read_file (the same regardless of
                             fraction membership, since every fraction reads the whole file).
    """
    counts = [0] * len(contexts)
    total_length = 0

    sys.stderr.write(f"Reading in {read_file}\n")
    for record in pyfastx.Fastq(read_file, build_index=False):
        name, seq, qual = record
        total_length += len(seq)
        trimmed_name = trim_read_id(name)

        for i, ctx in enumerate(contexts):
            taxon_id = ctx["read_map"].get(trimmed_name)
            if ctx["inverse"]:
                if taxon_id is not None:
                    for taxon in ctx["subtaxa_map"][taxon_id]:
                        update_summary_with_record(
                            taxon, record, ctx["out_counts"], ctx["read_stats"]
                        )
                    continue
                else:
                    counts[i] += 1
                    update_summary_with_record(
                        "other", record, ctx["out_counts"], ctx["read_stats"]
                    )
                    ctx["writer"].write("other", record)
            else:
                if taxon_id is None:
                    continue
                for taxon in ctx["subtaxa_map"][taxon_id]:
                    counts[i] += 1
                    update_summary_with_record(
                        taxon, record, ctx["out_counts"], ctx["read_stats"]
                    )
                    ctx["writer"].write(taxon, record)

    for ctx in contexts:
        ctx["writer"].close()

    return counts, total_length


def process_read_files_multi(
    fraction_specs: list,
    filetype: str,
    read_file_1: Path,
    read_file_2: Path = None,
) -> tuple[list, int]:
    """
    Like process_read_files, but for several independent fractions sharing a single pass over
    (paired/unpaired) read_files, instead of one pass per fraction.

    Args:
        fraction_specs (list): A list of dicts, one per fraction, each with keys 'prefixes'
                                (dict), 'read_map' (dict) and 'subtaxa_map' (dict) - the same
                                per-fraction inputs process_read_files takes for a single fraction.
        filetype (str): "fasta" or "fastq" depending on input filetype.
        read_file_1 (Path): Name of read file.
        read_file_2 (Path): Name of read file pair (if it exists).

    Returns:
        results (list): A list of (out_counts, read_stats, filenames) tuples, one per fraction,
                         in the same order as fraction_specs.
        total_length (int): Total length of all reads across the input file(s) - the same value
                             for every fraction, computed once here instead of once per fraction.
    """
    contexts = []
    for spec in fraction_specs:
        filenames = setup_outfiles(read_file_2, spec["prefixes"], filetype)
        contexts.append(
            {
                "read_map": spec["read_map"],
                "subtaxa_map": spec["subtaxa_map"],
                "inverse": spec["inverse"],
                "out_counts": defaultdict(int),
                "read_stats": defaultdict(new_taxon_stats),
                "filenames": filenames,
                "writer": TaxonWriter(filenames, 0),
            }
        )

    forward_counts, total_length = file_iterator_multi(read_file_1, contexts)

    if read_file_2:
        for ctx in contexts:
            ctx["writer"] = TaxonWriter(ctx["filenames"], 1)
        reverse_counts, reverse_length = file_iterator_multi(read_file_2, contexts)
        total_length += reverse_length
        for forward_count, reverse_count in zip(forward_counts, reverse_counts):
            if forward_count != reverse_count and (
                forward_count == 0 or reverse_count == 0
            ):
                sys.stderr.write(
                    "ERROR: No reads found for one of the file pair: extracted %i an %i reads respectively"
                    % (forward_count, reverse_count)
                )
                sys.exit(7)

    return (
        [
            (ctx["out_counts"], ctx["read_stats"], ctx["filenames"])
            for ctx in contexts
        ],
        total_length,
    )


def process_read_files(
    prefixes: str,
    filetype: str,
    read_map: dict,
    subtaxa_map: dict,
    read_file_1: Path,
    read_file_2: Path = None,
    inverse: bool = False,
) -> tuple[dict, dict, dict, dict]:
    """
    Iterate through (paired/unpaired) read_files and save the relevant reads, collecting summary statistics for these reads.

    Args:
        prefixes (dict): A dictionary from the key (usually taxon_id) to the prefix of the output file for that key.
        filetype (str): "fasta" or "fastq" depending on input filetype.
        read_map (dict): Dictionary from read ID to list of taxon ID associated with it.
        subtaxa_map (dict): Dictionary from taxon ID (output by read map) to list of taxon ID associated with
                            out files to be updated.
        read_file_1 (Path): Name of read file.
        read_file_2 (Path): Name of read file pair (if it exists).
        inverse (bool): If True, exclude all reads which match the taxon IDs in the read_map. If False, select them.

    Returns:
        out_counts (dict): Taxon ID to read count for that taxon.
        read_stats (dict): Taxon ID to running {"count", "sum_qual", "sum_len"} accumulator for that taxon.
        filenames (dict): Taxon ID to list of output filenames.
    """
    out_counts = defaultdict(int)
    read_stats = defaultdict(new_taxon_stats)

    filenames = setup_outfiles(read_file_2, prefixes, filetype)

    # Each writer is scoped to a single read file (forward or reverse) so its
    # LRU handle pool doesn't have to span both passes.
    writer_1 = TaxonWriter(filenames, 0)
    forward_count, total_length = file_iterator(
        read_file_1,
        read_map,
        subtaxa_map,
        inverse,
        out_counts,
        read_stats,
        writer_1,
    )

    reverse_count = 0
    if read_file_2:
        writer_2 = TaxonWriter(filenames, 1)
        reverse_count, reverse_length = file_iterator(
            read_file_2,
            read_map,
            subtaxa_map,
            inverse,
            out_counts,
            read_stats,
            writer_2,
        )
        total_length += reverse_length
        if forward_count != reverse_count and (
            forward_count == 0 or reverse_count == 0
        ):
            sys.stderr.write(
                "ERROR: No reads found for one of the file pair: extracted %i an %i reads respectively"
                % (forward_count, reverse_count)
            )
            sys.exit(7)
    return (out_counts, read_stats, filenames, total_length)


def generate_summary(
    lists_to_extract,
    entries,
    prefix,
    out_counts,
    read_stats,
    filenames,
    total_length,
    include_unclassified=False,
    short=False,
    write_total_length=True,
):
    """
    Generate a summary JSON file, including information about each taxon ID and corresponding read statistics.

    Args:
        lists_to_extract (list): A list of taxon ID to include in the summary.
        entries (dict): A dictionary containing information about the name/rank associated with each taxon ID, either
                        inferred from the Kraken report or NCBI taxonomy.
        prefix (str): Prefix for the summary file.
        out_counts (dict): Taxon ID to read count for that taxon.
        read_stats (dict): Taxon ID to running {"count", "sum_qual", "sum_len"} accumulator for that taxon.
        filenames (dict): Taxon ID to list of output filenames.
        include_unclassified (bool): True if each output file includes unclassified reads as well as those associated
                                     with the taxon ID.
        short (bool): Output short form summary, excluding taxon information.
        write_total_length (bool): Whether to write total_length.json. total_length is the length
                                    of the whole input file, the same value regardless of which
                                    fraction/report is being summarised - when generate_summary is
                                    called once per fraction sharing the same read file(s), only
                                    one of those calls needs to write it.
        filenames (dict): Taxon ID to list of output filenames.
    """
    sys.stderr.write("Write summary\n")
    summary = []
    for taxon_id in out_counts:
        if out_counts[taxon_id] == 0:
            sys.stderr.write(f"No reads extracted  for taxon_id {taxon_id}\n")
            continue

        taxon_stats = read_stats[taxon_id]
        avg_qual = (
            taxon_stats["sum_qual"] / taxon_stats["count"]
            if taxon_stats["count"]
            else 0
        )
        mean_len = (
            taxon_stats["sum_len"] / taxon_stats["count"]
            if taxon_stats["count"]
            else 0
        )
        qc_metrics = {
            "num_reads": out_counts[taxon_id],
            "avg_qual": avg_qual,
            "mean_len": mean_len,
            "total_len": taxon_stats["sum_len"],
        }

        if short:
            summary.append(
                {
                    "filenames": filenames[taxon_id],
                    "qc_metrics": qc_metrics,
                    "includes_unclassified": include_unclassified,
                }
            )
        else:
            summary.append(
                {
                    "human_readable": entries[taxon_id].name,
                    "taxon_id": taxon_id,
                    "tax_level": entries[taxon_id].simple_rank,
                    "filenames": filenames[taxon_id],
                    "qc_metrics": qc_metrics,
                    "includes_unclassified": include_unclassified,
                }
            )

    with open(f"{prefix}_summary.json", "w") as f:
        json.dump(summary, f)

    if write_total_length:
        with open("total_length.json", "w") as f:
            json.dump({"total_len": total_length}, f)

    for taxon_id in read_stats:
        taxon_stats = read_stats[taxon_id]
        if taxon_stats["count"] and taxon_stats["sum_qual"] / taxon_stats["count"] > 60:
            sys.stderr.write(
                f"Mean quality score for taxon_id {taxon_id} is "
                f"{taxon_stats['sum_qual'] / taxon_stats['count']} > 60. This seems "
                "unlikely! Please report this example to the DIPI group for investigation\n"
            )
            sys.exit(1)
