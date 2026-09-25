#!/usr/bin/env python3
"""
transfer_mod_tags.py

Copy MM/Mm and ML tags from a source BAM (demuxed basecalled BAM with mod tags)
to a mapped BAM (aligned, sorted), matching by read name (QNAME).

Usage:
  transfer_mod_tags.py --source demux/sample.bam --mapped mapped/sample.sorted.bam --out mapped/sample.sorted.mmtags.bam

Requires: pysam
"""
import argparse, sys
try:
    import pysam
except Exception as e:
    sys.stderr.write("ERROR: pysam is required. Install with pip install pysam or conda install -c bioconda pysam\n")
    raise

def has_mod_tags(bam_path, max_reads=10000):
    try:
        with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bamfile:
            for i, read in enumerate(bamfile.fetch(until_eof=True)):
                if read.has_tag("MM"):
                    return True
                if i > max_reads:
                    break
    except Exception as e:
        print(f"[WARN] Could not read BAM {bam_path}: {e}")
    return False


def build_tag_map(src_bam, tags_to_copy=("MM","ML","Mm")):
    tagmap = {}
    count = 0
    with pysam.AlignmentFile(src_bam, "rb", check_sq=False) as sf:
        for r in sf:
            tags = {}
            for tag_name in tags_to_copy:
                try:
                    val = r.get_tag(tag_name)
                except KeyError:
                    continue
                tags[tag_name] = val
            if tags:
                tagmap[r.query_name] = tags
                count += 1
    sys.stderr.write(f"[INFO] Collected modification tags for {count} reads from {src_bam}\n")
    return tagmap

def transfer_tags(mapped_bam, tagmap, out_bam):
    n_total = 0
    n_added = 0
    with pysam.AlignmentFile(mapped_bam, "rb") as mf, \
         pysam.AlignmentFile(out_bam, "wb", template=mf) as of:
        for r in mf:
            n_total += 1
            qname = r.query_name
            if qname in tagmap:
                for tname, tval in tagmap[qname].items():
                    # overwrite or set
                    r.set_tag(tname, tval, replace=True)
                n_added += 1
            of.write(r)
    sys.stderr.write(f"[INFO] Processed {n_total} alignments; added tags to {n_added} reads -> {out_bam}\n")
    return n_total, n_added

def index_bam(bamfile):
    pysam.index(bamfile)
    sys.stderr.write(f"[INFO] Indexed {bamfile}\n")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", required=True, help="Source demux BAM with mod tags (MM/ML)")
    p.add_argument("--mapped", required=True, help="Mapped sorted BAM (no mod tags)")
    p.add_argument("--out", required=True, help="Output BAM path (will be written, then indexed)")
    p.add_argument("--tags", nargs="+", default=["MM","ML","Mm"], help="Tags to copy (default: MM ML Mm)")
    args = p.parse_args()

    tagmap = build_tag_map(args.source, tags_to_copy=tuple(args.tags))
    transfer_tags(args.mapped, tagmap, args.out)
    index_bam(args.out)

if __name__ == "__main__":
    main()
