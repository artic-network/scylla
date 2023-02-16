#!/usr/bin/env python

import argparse
import sys
import mappy as mp
from collections import defaultdict


def load_manifest(log, path, preset):
    manifest = {
        "preset": preset,
        "references": [
        ],
    }
    manifest_fh = open(path)
    for line_i, line in enumerate(manifest_fh):
        fields = line.strip().split() # split on any whitespace if you have whitespace in your ref name you have bigger problems

        if line[0] == '#':
            continue

        if len(fields) < 3:
            sys.stderr.write("[FAIL] Manifest did not contain a third column mapping a reference to a preset\n")
            sys.stderr.write("       Consult the README to ensure you are using a manifest suitable for dehumaniser >= 0.9.0\n")
            sys.exit(78) # EX_CONFIG

        if fields[2] != preset:
            continue

        manifest["references"].append({
            "name": fields[0],
            "path": fields[1],
        })

    if len(manifest["references"]) == 0:
        sys.stderr.write("[FAIL] Manifest did not contain any references for preset=%s\n" % preset)
        sys.stderr.write("       Consult the README to ensure your manifest is correctly configured and for\n")
        sys.stderr.write("       instructions on how to build your own indexes if needed\n")
        sys.exit(65) # EX_DATAERR
    else:
        log.write("[NOTE] Detected %d references in manifest for preset=%s\n" % (len(manifest["references"]), preset))
    manifest_fh.close()
    return manifest


def dh_fastx_minimap(log, manifest, args):

    aligners = []
    for ref_i, ref_manifest in enumerate(manifest["references"]):
        aligners.append( mp.Aligner(ref_manifest["path"], preset=manifest["preset"]) )

    input_fastx = args.input

    if args.clean == "-":
        clean_fq = sys.stdout
        sys.stderr.write("[INFO] Writing to STDOUT\n")
    else:
        clean_fq = open(args.clean, 'w')
        sys.stderr.write("[INFO] Writing FASTX %s\n" % (args.clean))

    dirty_fq = None
    if args.dirty:
        dirty_fq = open(args.dirty, 'w')
        log.write("[INFO] Writing bad reads to FASTX %s\n" % (args.dirty))

    n_seqs = 0
    bad_seqs = 0
    for name, seq, qual in mp.fastx_read(input_fastx):
        n_seqs += 1
        if n_seqs % 10000 == 0:
            log.write("[INFO] Processed %d reads and dropped %d sequences\n" % (n_seqs ,bad_seqs))

        filter = False
        for ref_i, ref_manifest in enumerate(manifest["references"]):
            for hit in aligners[ref_i].map(seq):
                if dirty_fq:
                    if qual is None:
                        out_read = ">%s\n%s\n" % (name, seq)
                    else:
                        out_read = "@%s\n%s\n+\n%s\n" % (name,seq,qual)
                    dirty_fq.write(out_read)
                filter = True
                break

        if filter:
            bad_seqs += 1
            if args.flag_number and bad_seqs >= args.flag_number:
                sys.exit("[WARNING] Found %d bad sequences - quitting without finishing\n" % bad_seqs)
        else:
            if qual is None:
                out_read = ">%s\n%s\n" % (name, seq)
            else:
                out_read = "@%s\n%s\n+\n%s\n" % (name,seq,qual)
            clean_fq.write(out_read)
    clean_fq.close()
    if dirty_fq:
        dirty_fq.close()

    log.write("[INFO] Dropped %d sequences out of %d\n" % (bad_seqs, n_seqs))

def search_kraken_assignments_file(assignments_file, log, r_taxid="9606"):
    reads = []
    n_seqs = 0
    with open(assignments_file,"r") as f:
        for line in f:
            n_seqs += 1
            if n_seqs % 10000 == 0:
                log.write("[INFO] Processed %d reads and dropped %d sequences\n" % (n_seqs ,len(reads)))
            if r_taxid not in line:
                continue
            c_status, read, taxid, n, hits = line.strip().split("\t")
            key_count = defaultdict(int)
            for hit in hits.split():
                t,t_count = hit.split(":")
                key_count[t] += int(t_count)
            del key_count['0']
            if len(key_count) == 1:
                reads.append(read)
            else:
                bad, other = 0,0
                for k in key_count:
                    if k == r_taxid:
                        bad += key_count[k]
                    else:
                        other += key_count[k]
                if bad >= other:
                    reads.append(read)
    log.write("[INFO] Found %d sequences out of %d\n" % (len(reads), n_seqs))
    return reads

def scrub_by_read_id(in_fastx, out_fastx, bad_read_ids, log, bad_fastx=None):
    log.write("[INFO] Filter read file\n" )
    clean = open(out_fastx,"w")
    if bad_fastx:
        dirty = open(bad_fastx,"w")

    for name, seq, qual in mp.fastx_read(in_fastx):
        if name in bad_read_ids:
            if bad_fastx:
                if qual is None:
                    out_read = ">%s\n%s\n" % (name, seq)
                else:
                    out_read = "@%s\n%s\n+\n%s\n" % (name,seq,qual)
                dirty.write(out_read)
            continue
        elif qual is None:
            out_read = ">%s\n%s\n" % (name, seq)
        else:
            out_read = "@%s\n%s\n+\n%s\n" % (name,seq,qual)
        clean.write(out_read)
    clean.close()

def dh_fastx_kraken(log, args):
    read_ids = search_kraken_assignments_file(args.kraken_assignments, log)
    scrub_by_read_id(args.input, args.clean, read_ids, log, args.dirty)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="input dirty file")

    parser.add_argument("--manifest", help="reference manifest")
    parser.add_argument("--preset", help="mappy aligner preset")

    parser.add_argument("--kraken_assignments", help="kraken assignments file")

    parser.add_argument("-o", "--clean", help="output clean file [default -]", default="-")
    parser.add_argument("-d", "--dirty", help="output dirty file [default None]", default=None)

    parser.add_argument("-f", "--flag_number", type=int, help="warn and quit if more than this number of reads have been found [default None]", default=None)

    parser.add_argument("--log", help="log path [default <input>.dehumanizer.log.txt]", default=None)

    args = parser.parse_args()

    if not args.log:
        log = open(args.input + ".dehumanizer.log.txt", 'w')
    else:
        log = open(args.log, 'w')

    if args.manifest:
        if not args.preset:
            sys.stderr.write("Preset must be specified with --preset in this mode.\n")
            sys.exit()

        manifest = load_manifest(log, args.manifest, args.preset)

        dh_fastx_minimap(log, manifest, args)
    elif args.kraken_assignments:
        dh_fastx_kraken(log, args)
    else:
        sys.stderr.write("Must provide either --kraken or --manifest.\n")
        sys.exit()

    log.close()

if __name__ == "__main__":
    main()