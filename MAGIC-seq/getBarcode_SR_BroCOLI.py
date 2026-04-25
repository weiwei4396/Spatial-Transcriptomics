#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parser.add_argument("-i", required=True, help="input R1 fastq.gz")
parser.add_argument("-I", required=True, help="input R2 fastq.gz")
parser.add_argument("-x", default="", help="input barcodeX whitelist")
parser.add_argument("-y", default="", help="input barcodeY whitelist")
parser.add_argument("-z", default="", help="input barcodeZ whitelist")
parser.add_argument("-m", default=3, help="how many barcodes?")
parser.add_argument("-o", default="./", help="output folder path")
"""

import sys
import os
import time
import gzip
from itertools import zip_longest
import argparse


NUCS = 'ACGTN'


def smart_open(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def load_whitelist(path):
    barcodes = []
    with open(path) as f:
        for ln in f:
            ln = ln.rstrip('\n').rstrip('\r')
            if not ln:
                continue
            parts = ln.split('\t')
            if len(parts) < 2:
                continue
            bc = parts[1].strip()
            if bc:
                barcodes.append(bc)
    return barcodes


def build_lookup(whitelist):
    """
    返回 {query_barcode: corrected_barcode}。
    - 精确匹配: bc -> bc
    - 1 次替换邻居: variant -> bc
    - 模糊邻居 (同一个 variant 指向多个不同 bc) 会被剔除, 避免错误纠正
    """
    wl = set(whitelist)
    lookup = {bc: bc for bc in wl}
    ambiguous = set()
    for bc in wl:
        for i in range(len(bc)):
            prefix, orig, suffix = bc[:i], bc[i], bc[i+1:]
            for n in NUCS:
                if n == orig:
                    continue
                variant = prefix + n + suffix
                if variant in wl:
                    # variant 自己就是白名单里的条码，它的精确匹配优先
                    continue
                if variant in ambiguous:
                    continue
                prev = lookup.get(variant)
                if prev is None:
                    lookup[variant] = bc
                elif prev != bc:
                    # 两个白名单条码各自编辑距离 1 都能到 variant，歧义 → 丢弃
                    del lookup[variant]
                    ambiguous.add(variant)
    return lookup


def grouper(iterable, n):
    """每 n 行打包成一组"""
    args = [iter(iterable)] * n
    return zip_longest(*args)


def process(fq1, fq2, out1, out2, lookup_x, lookup_y,
            sx, ex, sy, ey, su, eu, min_qual=10, has_qual=True, lookup_z=None, sz=None, ez=None,
            write_discarded=False, discard1=None, discard2=None):
    stats = dict(total=0, mismatch_x=0, mismatch_y=0, mismatch_z=0, umi_del=0, kept=0)
    qual_cutoff = min_qual + 33  # 直接比 ASCII，省掉 ord() -33
    has_z = lookup_z is not None
    countTimes = 0

    d1 = open(discard1, 'w') if write_discarded else None
    d2 = open(discard2, 'w') if write_discarded else None
    def write_bad_reads(rec1, rec2):
        d1.writelines(rec1)
        d2.writelines(rec2)

    try:
        with smart_open(fq1) as f1, smart_open(fq2) as f2, \
            open(out1, 'w') as o1, open(out2, 'w') as o2:
            
            for rec1, rec2 in zip(grouper(f1, 4), grouper(f2, 4)):
                countTimes += 1
                if countTimes % 1000000 == 0:
                    print(f"have processed {countTimes/1000000}M reads...")
                if rec1[0] is None:
                    break
                stats['total'] += 1

                seq1 = rec1[1].rstrip('\n')
                qual1 = rec1[3].rstrip('\n')

                # 1) UMI 质量
                if has_qual:
                    umi_q = qual1[su:eu]
                    if any(ord(c) < qual_cutoff for c in umi_q):
                        stats['umi_del'] += 1
                        if write_discarded:
                            write_bad_reads(rec1, rec2)
                        continue

                # 2) X 条码
                bx = seq1[sx:ex]
                cx = lookup_x.get(bx)
                if cx is None:
                    stats['mismatch_x'] += 1
                    if write_discarded:
                        write_bad_reads(rec1, rec2)
                    continue

                # 3) Y 条码
                by = seq1[sy:ey]
                cy = lookup_y.get(by)
                if cy is None:
                    stats['mismatch_y'] += 1
                    if write_discarded:
                        write_bad_reads(rec1, rec2)
                    continue

                # 4) Z 条码（可选）
                if has_z:
                    bz = seq1[sz:ez]
                    cz = lookup_z.get(bz)
                    if cz is None:
                        stats['mismatch_z'] += 1
                        if write_discarded:
                            write_bad_reads(rec1, rec2)
                        continue

                # # 5) 统一修正序列
                # modified = False
                # if cx != bx:
                #     seq1 = seq1[:sx] + cx + seq1[ex:]
                #     modified = True
                # if cy != by:
                #     seq1 = seq1[:sy] + cy + seq1[ey:]
                #     modified = True
                # if has_z and cz != bz:
                #     seq1 = seq1[:sz] + cz + seq1[ez:]
                #     modified = True
                # if modified:
                #     rec1 = (rec1[0], seq1 + '\n', rec1[2], rec1[3])

                # 5) 重构 read1：barcode_x + barcode_y + (barcode_z) + umi
                umi_seq = seq1[su:eu]
                if has_z:
                    new_seq1  = cx + cy + cz + umi_seq
                    new_qual1 = qual1[sx:ex] + qual1[sy:ey] + qual1[sz:ez] + qual1[su:eu]
                else:
                    new_seq1  = cx + cy + umi_seq
                    new_qual1 = qual1[sx:ex] + qual1[sy:ey] + qual1[su:eu]

                rec1 = (rec1[0], new_seq1 + '\n', rec1[2], new_qual1 + '\n')

                o1.writelines(rec1)
                o2.writelines(rec2)
                stats['kept'] += 1
    finally:
        if d1 is not None:
            d1.close()
        if d2 is not None:
            d2.close()
    return stats

def strip_fq_ext(path):
    """去掉目录，再去掉 .fq / .fastq / .fq.gz / .fastq.gz 后缀"""
    name = os.path.basename(path)
    for ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
        if name.endswith(ext):
            return name[:-len(ext)]
    return name.rsplit('.', 1)[0]

def main():

    parser = argparse.ArgumentParser(description="get reads barcode!")

    parser.add_argument("-i", required=True, help="input R1 fastq.gz")
    parser.add_argument("-I", required=True, help="input R2 fastq.gz")
    parser.add_argument("-x", default="", help="input barcodeX whitelist")
    parser.add_argument("-y", default="", help="input barcodeY whitelist")
    parser.add_argument("-z", default="", help="input barcodeZ whitelist")
    parser.add_argument("-m", default=3, help="how many barcodes?")
    parser.add_argument("-o", default="./", help="output folder path")
    parser.add_argument("--write_discarded", action="store_true",
                        help="output discarded read pairs to fastq files")

    args = parser.parse_args()

    fq1, fq2 = args.i, args.I
    mode = int(args.m)
    outputPath = args.o
    write_discarded = args.write_discarded

    if mode == 2:
        wl_x_path, wl_y_path = args.x, args.y
        sx, ex = 0, 8
        sy, ey = 26, 34
        su, eu = 34, 46

    elif mode == 3:
        wl_x_path, wl_y_path, wl_z_path = args.x, args.y, args.z
        lookup_z = build_lookup(load_whitelist(wl_z_path))
        sx, ex = 0, 8
        sy, ey = 26, 34
        sz, ez = 52, 60
        su, eu = 60, 72
    
    t0 = time.time()
    print("--- *** Now we start to count *** ---")
    lookup_x = build_lookup(load_whitelist(wl_x_path))
    lookup_y = build_lookup(load_whitelist(wl_y_path))
    t_idx = time.time() - t0


    prefix1 = strip_fq_ext(fq1)
    prefix2 = strip_fq_ext(fq2)
    out1 = os.path.join(outputPath, prefix1 + "_trim.fq")
    out2 = os.path.join(outputPath, prefix2 + "_trim.fq")
    discard1 = os.path.join(outputPath, prefix1 + "_discarded.fq")
    discard2 = os.path.join(outputPath, prefix2 + "_discarded.fq")

    if mode == 2:
        stats = process(fq1, fq2, out1, out2, lookup_x, lookup_y,
                        sx, ex, sy, ey, su, eu,
                        write_discarded=write_discarded, discard1=discard1, discard2=discard2)
    elif mode == 3:
        stats = process(fq1, fq2, out1, out2, lookup_x, lookup_y,
                        sx, ex, sy, ey, su, eu,
                        lookup_z=lookup_z, sz=sz, ez=ez,
                        write_discarded=write_discarded, discard1=discard1, discard2=discard2)  
      
    run_min = (time.time() - t0) / 60

    print(f"before trim: {stats['total']}")
    print(f"after trim: {stats['kept']}, ratio: {stats['kept']/stats['total']*100}%")
    print(f"reads discarded in step1--mismatch_X: {stats['mismatch_x']}, ratio: {stats['mismatch_x']/stats['total']*100}%")
    print(f"reads discarded in step1--mismatch_Y: {stats['mismatch_y']}, ratio: {stats['mismatch_y']/stats['total']*100}%")
    if mode == 3:
        print(f"reads discarded in step1--mismatch_Z: {stats['mismatch_z']}, ratio: {stats['mismatch_z']/stats['total']*100}%")
    print(f"reads discarded in step2--low quality umi: {stats['umi_del']}, ratio: {stats['umi_del']/stats['total']*100}%")
    print(f"index build time: {t_idx:.2f} s\n")
    print(f"qc control completed, run time: {run_min:.2f} min\n")


    # with open(prefix1 + '_trim_report.txt', 'w') as r:
    #     r.write(f"before trim: {stats['total']}\n")
    #     r.write(f"after trim: {stats['kept']}\n")
    #     r.write(f"reads discarded in step1--mismatch_X: {stats['mismatch_x']}\n")
    #     r.write(f"reads discarded in step1--mismatch_Y: {stats['mismatch_y']}\n")
    #     r.write(f"reads discarded in step2--low quality umi: {stats['umi_del']}\n")
    #     r.write(f"index build time: {t_idx:.2f} s\n")
    #     r.write(f"qc control completed, run time: {run_min:.2f} min\n")

if __name__ == '__main__':
    main()


























