#!/usr/bin/env python
# coding=utf-8

"""
BLAST 污染检测脚本 v2
——基于 BLAST+ 内置 taxdb 直接输出 taxid / 物种名，
  不再依赖 accession2taxid、names.dmp 这类外部映射表。

用法:
    # 已有 fa 文件:
    python 1_blast_species_v2.py -f sample.fa -d /path/to/core_nt

    # 从 fastq 开始 (自动抽前 5000 条):
    python 1_blast_species_v2.py -f sample.fa -q sample.fastq -d /path/to/core_nt
"""
import os
import re
import sys
import time
import argparse
import subprocess
from collections import Counter


# ---------- 工具函数 ----------
def count_fa_reads(fa):
    """统计 fa 文件里 read 的总条数（> 开头的行数）。"""
    n = 0
    with open(fa, "r") as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def fastq_to_fa(fastq, fa, n_reads=5000):
    """从 fastq 头部抽 n_reads 条 read 转成 fasta。
       自动处理 .gz 压缩格式；用 awk 按行号判定（比 sed 匹配 @ 更可靠）。"""
    n_lines = n_reads * 4

    # 根据后缀自动选择读取方式: .gz 用 zcat, 明文用 cat
    if fastq.endswith(".gz"):
        reader = "zcat {fq}".format(fq=fastq)
    else:
        reader = "cat {fq}".format(fq=fastq)

    cmd = (
        "{reader} | head -n {nl} | "
        "awk '{{ if (NR%4==1){{ sub(/^@/, \">\"); print $1 }} "
        "else if (NR%4==2){{ print $0 }} }}' > {fa}"
    ).format(reader=reader, nl=n_lines, fa=fa)
    subprocess.run(cmd, shell=True, check=True)


def ensure_taxdb(db_path):
    """检查 taxdb.btd/.bti 是否在数据库目录下，
       没有的话给出清晰的提示（光有 core_nt.* 是不够的）。"""
    db_dir = os.path.dirname(db_path) or "."
    needed = ["taxdb.btd", "taxdb.bti"]
    missing = [f for f in needed if not os.path.exists(os.path.join(db_dir, f))]
    if missing:
        sys.stderr.write(
            "WARNING: {dir} 下缺少 {miss}。\n"
            "         没有这两个文件, BLAST 的 sscinames 字段会全部是 'N/A'。\n"
            "         到 https://ftp.ncbi.nlm.nih.gov/blast/db/ 下载 taxdb.tar.gz "
            "解压到该目录即可。\n".format(dir=db_dir, miss=missing)
        )
    # 让 BLAST 能定位到 taxdb
    os.environ["BLASTDB"] = db_dir + os.pathsep + os.environ.get("BLASTDB", "")


def run_blast(fa, blast_tsv, db, threads=20, evalue=1e-5):
    """运行 blastn, 直接输出带 taxid 和 sscinames 的 tsv。"""
    outfmt = "6 qseqid sseqid staxids sscinames scomnames pident length evalue bitscore"
    cmd = [
        "blastn",
        "-query", fa,
        "-out", blast_tsv,
        "-db", db,
        "-max_target_seqs", "1",
        "-outfmt", outfmt,
        "-num_threads", str(threads),
        "-evalue", str(evalue),
    ]
    print("Running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_and_summarize(blast_tsv, total_reads, out_path):
    """解析 tsv → 每条 read 保留最佳 hit → 按物种统计占比。"""
    cols = ["qseqid", "sseqid", "staxids", "sscinames", "scomnames",
            "pident", "length", "evalue", "bitscore"]

    # 每条 read 只保留 evalue 最小 / bitscore 最大的那行
    # (max_target_seqs=1 通常已经够了，这里再兜一层)
    best = {}
    with open(blast_tsv, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(cols):
                continue
            row = dict(zip(cols, parts))
            qid = row["qseqid"]
            ev = float(row["evalue"])
            bs = float(row["bitscore"])
            if qid not in best:
                best[qid] = (ev, bs, row)
            else:
                old_ev, old_bs, _ = best[qid]
                if ev < old_ev or (ev == old_ev and bs > old_bs):
                    best[qid] = (ev, bs, row)

    # 统计物种。sscinames / staxids 可能是 "A;B" 形式，取第一个
    species_list = []
    for _, _, row in best.values():
        name = row["sscinames"].strip()
        if not name or name == "N/A":
            name = "Unknown"
        name = name.split(";")[0].strip()
        species_list.append(name)

    hit_reads = len(species_list)
    counts = Counter(species_list)
    count_list = sorted(counts.items(), key=lambda x: x[1], reverse=True)

    with open(out_path, "w") as out:
        out.write("Name\tHit_reads\tpercentage_total\tpercentage_hits\n")
        for name, n in count_list:
            p1 = (100.0 * n / total_reads) if total_reads else 0.0
            p2 = (100.0 * n / hit_reads)   if hit_reads   else 0.0
            out.write("{name}\t{n}\t{p1:.2f}%\t{p2:.2f}%\n".format(
                name=name, n=n, p1=p1, p2=p2))

    return hit_reads


# ---------- 主流程 ----------

def main():
    parser = argparse.ArgumentParser(
        description="BLAST 污染检测 (新版: 基于 BLAST+ 内置 taxdb)"
    )
    parser.add_argument("-f", "--fa", default=None,
                        help="输入/输出 fa 文件路径。若只提供 --fastq 而不给本项,"
                             " 会自动从 fastq 文件名派生 (sample.fastq -> sample.fa)")
    parser.add_argument("-q", "--fastq", default=None,
                        help="可选: 原始 fastq, 会抽前 N 条 read 转为 fa")
    parser.add_argument("-d", "--db", required=True,
                        help="BLAST 数据库路径前缀, 如 /data/.../core/core_nt")
    parser.add_argument("-t", "--threads", type=int, default=20,
                        help="BLAST 线程数 (默认 20)")
    parser.add_argument("-e", "--evalue", type=float, default=1e-5,
                        help="BLAST evalue 阈值 (默认 1e-5)")
    parser.add_argument("-n", "--n-reads", type=int, default=5000,
                        help="从 fastq 抽取的 read 条数 (默认 5000)")
    args = parser.parse_args()

    # 参数校验: 至少要提供 -f 或 -q 其中一个
    if not args.fa and not args.fastq:
        parser.error("必须至少提供 -f/--fa 或 -q/--fastq 中的一个")

    # 若只提供了 fastq, 自动从 fastq 文件名派生 fa 路径
    # 支持 .fastq / .fq / .fastq.gz / .fq.gz 这几种常见后缀
    if args.fa:
        fa = args.fa
    else:
        fa = re.sub(r"\.(fastq|fq)(\.gz)?$", "", args.fastq) + ".fa"
        print("未指定 -f, 自动派生 fa 路径: {f}".format(f=fa))

    prefix = re.sub(r"\.fa$", "", fa)
    blast_tsv = prefix + "_blast_result.tsv"
    final_out = prefix + "_blast_final_result.txt"

    # 1) 已有最终结果则跳过
    if os.path.exists(final_out):
        print("Pollution already done, passed pollution\n")
        return

    t0 = time.time()

    # 2) fastq → fa (可选)
    if args.fastq and not os.path.exists(fa):
        print("Converting fastq -> fa (first {n} reads)".format(n=args.n_reads))
        fastq_to_fa(args.fastq, fa, args.n_reads)

    if not os.path.exists(fa):
        sys.stderr.write("ERROR: fa file not found: {f}\n".format(f=fa))
        sys.exit(1)

    # 3) 检查/暴露 taxdb 路径
    ensure_taxdb(args.db)

    # 4) 数一下 reads 总数
    total_reads = count_fa_reads(fa)
    print("Total reads in fa: {n}".format(n=total_reads))

    # 5) 运行 BLAST (已有结果则跳过, 方便调试时重复解析)
    print("1_pollution analysis start!")
    if not os.path.exists(blast_tsv):
        run_blast(fa, blast_tsv, args.db,
                  threads=args.threads, evalue=args.evalue)
    else:
        print("BLAST output already exists, skip blastn: {t}".format(t=blast_tsv))

    # 6) 解析并汇总
    hit_reads = parse_and_summarize(blast_tsv, total_reads, final_out)
    print("Total reads: {t}, hit reads: {h}".format(t=total_reads, h=hit_reads))
    print("Result written to: {o}".format(o=final_out))

    elapsed = (time.time() - t0) / 60.0
    print(" pollution completed, run time: {:.2f} min\n".format(elapsed))


if __name__ == "__main__":
    main()






