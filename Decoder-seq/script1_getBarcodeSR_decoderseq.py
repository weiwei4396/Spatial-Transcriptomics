#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parser.add_argument("-i", required=True, help="input R1 fastq.gz")
parser.add_argument("-I", required=True, help="input R2 fastq.gz")
parser.add_argument("-x", default="", help="input barcode whitelist")
parser.add_argument("--bcx", default="", help="input barcodeX whitelist")
parser.add_argument("--bcy", default="", help="input barcodeY whitelist")
parser.add_argument("-o", default="./", help="output folder path")
parser.add_argument("--write_discarded", action="store_true",
                    help="output discarded read pairs to fastq files")
parser.add_argument("--with_umi", action="store_true",
                    help="umi filter")
"""

import os
import time
import gzip
import argparse
from itertools import zip_longest

NUCS = 'ACGTN'

def smart_open(path, mode='rt'):
    """Open plain text or gzipped file."""
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def load_decoderseq_whitelist(path, barcode_mode=1):
    """
    Load Decoder-seq whitelist.

    barcode_mode = 1: use column index 4 as Barcode X
    barcode_mode = 2: use column index 5 as Barcode Y
    """
    barcodes = []
    line_number = 0

    with open(path) as f:
        for ln in f:
            line_number += 1
            ln = ln.rstrip("\n").rstrip("\r")

            if not ln:
                continue

            parts = ln.split("\t")

            if barcode_mode == 1:
                if len(parts) < 5:
                    continue
                bc = parts[4].strip()
            elif barcode_mode == 2:
                if len(parts) < 6:
                    continue
                bc = parts[5].strip()
            else:
                raise ValueError("barcode_mode must be 1 or 2")

            if bc:
                barcodes.append(bc)

    if barcode_mode == 1:
        print(f"Barcode X whitelist loaded, total {len(barcodes)} barcodes.")
    elif barcode_mode == 2:
        print(f"Barcode Y whitelist loaded, total {len(barcodes)} barcodes.")

    return barcodes

def load_whitelist(path):
    barcodes = []
    with open(path) as f:
        for ln in f:
            ln = ln.rstrip('\n').rstrip('\r')
            if not ln:
                continue
            parts = ln.split('\t')
            bc = parts[0].strip()
            if bc:
                barcodes.append(bc)
    return barcodes

def build_lookup(whitelist):
    """
    Build barcode correction lookup.

    Return:
        {query_barcode: corrected_barcode}

    Includes:
        exact match:
            bc -> bc

        1-substitution correction:
            variant -> bc

    Ambiguous variants are removed.
    For example, if one observed variant can be corrected to two different
    whitelist barcodes, it will not be used.
    """
    wl = set(whitelist)
    lookup = {bc: bc for bc in wl}
    ambiguous = set()

    for bc in wl:
        for i in range(len(bc)):
            prefix = bc[:i]
            orig = bc[i]
            suffix = bc[i + 1:]

            for n in NUCS:
                if n == orig:
                    continue

                variant = prefix + n + suffix

                if variant in wl:
                    # Exact whitelist barcode has priority.
                    continue

                if variant in ambiguous:
                    continue

                prev = lookup.get(variant)

                if prev is None:
                    lookup[variant] = bc
                elif prev != bc:
                    # Ambiguous correction. Remove it.
                    if variant in lookup:
                        del lookup[variant]
                    ambiguous.add(variant)

    return lookup


def grouper(iterable, n):
    """Group FASTQ lines into records of n lines."""
    args = [iter(iterable)] * n
    return zip_longest(*args)


def write_fastq_record(handle, record):
    """Write one FASTQ record."""
    handle.writelines(record)


def process(fq1, fq2, out1, out2, lookup_x, lookup_y,
            sx, ex, sy, ey, su, eu, min_qual=10, has_qual=False, sz=None, ez=None,
            write_discarded=False, discard1=None, discard2=None):
    stats = dict(
        total=0,
        mismatch_x=0,
        mismatch_y=0,
        umi_del=0,
        kept=0,
        corrected_x=0,
        corrected_y=0,
        corrected_any=0,
        corrected_both_xy=0,
    )
    qual_cutoff = min_qual + 33  # 直接比 ASCII，省掉 ord() -33
    count_times = 0

    d1 = open(discard1, 'w') if write_discarded else None
    d2 = open(discard2, 'w') if write_discarded else None
    def write_bad_reads(rec1, rec2):
        if d1 is not None and d2 is not None:
            d1.writelines(rec1)
            d2.writelines(rec2)

    try:
        with smart_open(fq1) as f1, smart_open(fq2) as f2, \
                open(out1, "w") as o1, open(out2, "w") as o2:

            for rec1, rec2 in zip(grouper(f1, 4), grouper(f2, 4)):
                count_times += 1

                if count_times % 1000000 == 0:
                    print(f"have processed {count_times / 1000000:.0f}M reads...")

                if rec1[0] is None or rec2[0] is None:
                    break

                stats["total"] += 1

                seq1 = rec1[1].rstrip("\n")
                qual1 = rec1[3].rstrip("\n")

                # 1. Barcode X
                bx = seq1[sx:ex]
                cx = lookup_x.get(bx)

                if cx is None:
                    stats["mismatch_x"] += 1
                    if write_discarded:
                        write_bad_reads(rec1, rec2)
                    continue

                x_corrected = cx != bx

                # 2. Barcode Y
                by = seq1[sy:ey]
                cy = lookup_y.get(by)

                if cy is None:
                    stats["mismatch_y"] += 1
                    if write_discarded:
                        write_bad_reads(rec1, rec2)
                    continue

                y_corrected = cy != by

                # 3.Optional UMI quality filtering.
                if has_qual:
                    umi_q = qual1[su:eu]
                    if any(ord(c) < qual_cutoff for c in umi_q):
                        stats["umi_del"] += 1
                        if write_discarded:
                            write_bad_reads(rec1, rec2)
                        continue
                        

                # 4. Count corrected barcodes.
                if x_corrected:
                    stats["corrected_x"] += 1

                if y_corrected:
                    stats["corrected_y"] += 1

                if x_corrected or y_corrected:
                    stats["corrected_any"] += 1

                if x_corrected and y_corrected:
                    stats["corrected_both_xy"] += 1

                # 5. Reconstruct read1:
                # barcode_x + barcode_y + optional barcode_z + UMI
                umi_seq = seq1[su:eu]

                new_seq1 = cx + cy + umi_seq
                new_qual1 = qual1[sx:ex] + qual1[sy:ey] + qual1[su:eu]

                new_rec1 = (
                    rec1[0],
                    new_seq1 + "\n",
                    rec1[2],
                    new_qual1 + "\n",
                )

                o1.writelines(new_rec1)
                o2.writelines(rec2)
                stats["kept"] += 1
    
    finally:
        if d1 is not None:
            d1.close()
        if d2 is not None:
            d2.close()

    return stats

def strip_fq_ext(path):
    """Remove FASTQ suffix."""
    name = os.path.basename(path)

    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(ext):
            return name[:-len(ext)]

    return name.rsplit(".", 1)[0]

def pct(num, den):
    """Safe percentage."""
    if den == 0:
        return 0.0
    return num / den * 100

def main():

    parser = argparse.ArgumentParser(description="get reads barcode!")

    parser.add_argument("-i", required=True, help="input R1 fastq.gz")
    parser.add_argument("-I", required=True, help="input R2 fastq.gz")
    parser.add_argument("-x", default="", help="input barcodeX whitelist")
    parser.add_argument("--bcx", default="", help="input barcodeX whitelist")
    parser.add_argument("--bcy", default="", help="input barcodeX whitelist")
    parser.add_argument("-o", default="./", help="output folder path")
    parser.add_argument("--write_discarded", action="store_true",
                        help="output discarded read pairs to fastq files")
    parser.add_argument("--with_umi", action="store_true",
                        help="umi filter")

    args = parser.parse_args()

    fq1, fq2 = args.i, args.I
    outputPath = args.o
    write_discarded = args.write_discarded
    with_umi = args.with_umi

    wl_x_path = args.x
    sx, ex = 0, 8
    sy, ey = 38, 46
    su, eu = 46, 58
    
    t0 = time.time()
    print("--- *** Now we start to count *** ---")
    if args.bcx and args.bcy:
        lookup_x = build_lookup(load_whitelist(args.bcx))
        lookup_y = build_lookup(load_whitelist(args.bcy))
    else:
        lookup_x = build_lookup(load_decoderseq_whitelist(wl_x_path, barcode_mode=1))
        lookup_y = build_lookup(load_decoderseq_whitelist(wl_x_path, barcode_mode=2))
    t_idx = time.time() - t0

    prefix1 = strip_fq_ext(fq1)
    prefix2 = strip_fq_ext(fq2)

    out1 = os.path.join(outputPath, prefix1 + "_trim.fastq")
    out2 = os.path.join(outputPath, prefix2 + "_trim.fastq")
    discard1 = os.path.join(outputPath, prefix1 + "_discarded.fastq")
    discard2 = os.path.join(outputPath, prefix2 + "_discarded.fastq")

    stats = process(fq1, fq2, out1, out2, lookup_x, lookup_y,
                    sx, ex, sy, ey, su, eu, has_qual=with_umi,
                    write_discarded=write_discarded, discard1=discard1, discard2=discard2)

    run_min = (time.time() - t0) / 60

    print("\n========== Summary ==========")
    total = stats["total"]
    kept = stats["kept"]

    print(f"before trim: {total}")
    print(f"after trim: {kept}, ratio among total: {pct(kept, total):.4f}%")

    print(
        f"reads discarded in step1--mismatch_X: {stats['mismatch_x']}, "
        f"ratio among total: {pct(stats['mismatch_x'], total):.4f}%"
    )

    print(
        f"reads discarded in step1--mismatch_Y: {stats['mismatch_y']}, "
        f"ratio among total: {pct(stats['mismatch_y'], total):.4f}%"
    )

    print(
        f"reads discarded in step2--low quality UMI: {stats['umi_del']}, "
        f"ratio among total: {pct(stats['umi_del'], total):.4f}%"
    )

    print("\n========== Barcode correction ==========")

    print(
        f"barcode corrected X: {stats['corrected_x']}, "
        f"ratio among total: {pct(stats['corrected_x'], total):.4f}%, "
        f"ratio among kept: {pct(stats['corrected_x'], kept):.4f}%"
    )

    print(
        f"barcode corrected Y: {stats['corrected_y']}, "
        f"ratio among total: {pct(stats['corrected_y'], total):.4f}%, "
        f"ratio among kept: {pct(stats['corrected_y'], kept):.4f}%"
    )

    print(
        f"barcode corrected any X/Y: {stats['corrected_any']}, "
        f"ratio among total: {pct(stats['corrected_any'], total):.4f}%, "
        f"ratio among kept: {pct(stats['corrected_any'], kept):.4f}%"
    )

    print(
        f"barcode corrected both X and Y: {stats['corrected_both_xy']}, "
        f"ratio among total: {pct(stats['corrected_both_xy'], total):.4f}%, "
        f"ratio among kept: {pct(stats['corrected_both_xy'], kept):.4f}%"
    )

    print("\n========== Runtime ==========")

    print(f"index build time: {t_idx:.2f} s")
    print(f"qc control completed, run time: {run_min:.2f} min")

    print("\n========== Output ==========")
    print(f"output R1: {out1}")
    print(f"output R2: {out2}")

    if write_discarded:
        print(f"discarded R1: {discard1}")
        print(f"discarded R2: {discard2}")


if __name__ == '__main__':
    main()


























