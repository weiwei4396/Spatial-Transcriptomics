import argparse
import random
import pandas as pd
import matplotlib.pyplot as plt
import pysam
import numpy as np


def get_tag_or_none(read, tag):
    try:
        return read.get_tag(tag)
    except KeyError:
        return None


def calculate_saturation_from_bam(
    bam_path, out_tsv, out_png,
    fractions=None, seed=123, use_gene=True, min_mapq=0
):
    """
    Calculate sequencing saturation and median genes per spot from STARsolo BAM.
    Sequencing saturation = 1 - unique_molecules / valid_reads

    For each subsampling fraction:
        x axis = mean reads per spot
        y1 = sequencing saturation
        y2 = median genes per spot
    """

    if fractions is None:
        fractions = [i / 10 for i in range(1, 11)]

    fractions = sorted(fractions)

    # 初始化每个fraction下的read数量和molecule集合
    valid_reads = {f: 0 for f in fractions}
    molecules = {f: set() for f in fractions}

    # 每个fraction下, 每个spot检测到的gene集合;
    spot_genes = {f: {} for f in fractions}

    random.seed(seed)

    bam = pysam.AlignmentFile(bam_path, "rb")

    # n_total:扫描过的BAM记录总数; n_used:通过过滤后真正用于计算的read数;
    # all_spots:所有出现过的有效spot barcode集合;
    n_total = 0
    n_used = 0
    all_spots = set()

    # 遍历BAM文件中的每一条alignment record;
    # 表示即使BAM没有index, 也可以从头读到尾;
    for read in bam.fetch(until_eof=True):
        n_total += 1

        if read.is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        # cell barcode 和 umi barcode;
        cb = get_tag_or_none(read, "CB")
        ub = get_tag_or_none(read, "UB")

        if cb is None or ub is None:
            continue

        gx = None
        if use_gene:
            gx = get_tag_or_none(read, "GX")
            if gx is None:
                gx = get_tag_or_none(read, "GN")
            if gx is None:
                continue
            molecule = (cb, ub, gx)
        else:
            molecule = (cb, ub)

        n_used += 1
        all_spots.add(cb)

        r = random.random()

        for f in fractions:
            if r <= f:
                valid_reads[f] += 1
                molecules[f].add(molecule)

                # 只有 use_gene=True 时，才能统计 genes per spot
                if use_gene:
                    if cb not in spot_genes[f]:
                        spot_genes[f][cb] = set()
                    spot_genes[f][cb].add(gx)

    bam.close()

    # 至少有一条有效read的barcode;
    n_spots = len(all_spots)
    records = []

    for f in fractions:
        vr = valid_reads[f]
        um = len(molecules[f])

        if vr > 0:
            sat = 1 - um / vr
        else:
            sat = float("nan")

        mean_reads_per_spot = vr / n_spots if n_spots > 0 else float("nan")

        # 把所有有效spot都纳入计算, 没有检测到gene的spot记为0;
        if use_gene and n_spots > 0:
            genes_per_spot = [
                len(spot_genes[f].get(cb, set()))
                for cb in all_spots
            ]
            median_genes_per_spot = float(np.median(genes_per_spot))
        else:
            median_genes_per_spot = float("nan")

        records.append({
            "fraction": f,
            "valid_reads": vr,
            "unique_molecules": um,
            "sequencing_saturation": sat,
            "sequencing_saturation_percent": sat * 100,
            "mean_reads_per_spot": mean_reads_per_spot,
            "median_genes_per_spot": median_genes_per_spot
        })
        print(f"fraction:{f}")
        print(f"valid_reads:{vr}")
        print(f"unique_molecules:{um}")
        print(f"sequencing_saturation:{sat}")
        print(f"sequencing_saturation_percent:{sat * 100}%")
        print(f"mean_reads_per_spot:{mean_reads_per_spot}")
        print(f"median_genes_per_spot:{median_genes_per_spot}")

    df = pd.DataFrame(records)
    df.to_csv(out_tsv, sep="\t", index=False)

    # 画两个图
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    axes[0].plot(
        df["mean_reads_per_spot"],
        df["sequencing_saturation"],
        marker="o"
    )
    axes[0].set_xlabel("Mean Reads per Spot")
    axes[0].set_ylabel("Sequencing Saturation")
    axes[0].set_title("Sequencing Saturation")

    axes[1].plot(
        df["mean_reads_per_spot"],
        df["median_genes_per_spot"],
        marker="o"
    )
    axes[1].set_xlabel("Mean Reads per Spot")
    axes[1].set_ylabel("Median Genes per Spot")
    axes[1].set_title("Median Genes per Spot")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"Total BAM records scanned: {n_total}")
    print(f"Valid reads used: {n_used}")
    print(f"Number of spots: {n_spots}")
    print(f"Saved table to: {out_tsv}")
    print(f"Saved plot to: {out_png}")


# Sequencing saturation = 1 - unique_molecules / valid_reads;
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bam",
        required=True,
        help="STARsolo BAM file with CB/UB/GX tags"
    )
    parser.add_argument(
        "--out_tsv",
        required=True,
        help="Output TSV file"
    )
    parser.add_argument(
        "--out_png",
        required=True,
        help="Output PNG file"
    )
    parser.add_argument(
        "--fractions",
        default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0",
        help="Comma-separated read fractions"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed"
    )
    parser.add_argument(
        "--no_gene",
        action="store_true",
        help="Use CB+UB only instead of CB+UB+GX/GN"
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=0,
        help="Minimum mapping quality"
    )

    args = parser.parse_args()

    fractions = [float(x) for x in args.fractions.split(",")]
    if any(f <= 0 or f > 1 for f in fractions):
        raise ValueError("All fractions must be > 0 and <= 1.")

    calculate_saturation_from_bam(
        bam_path=args.bam,
        out_tsv=args.out_tsv,
        out_png=args.out_png,
        fractions=fractions,
        seed=args.seed,
        use_gene=not args.no_gene,
        min_mapq=args.min_mapq
    )


if __name__ == "__main__":
    main()


# python plot_sequencing_saturation_from_bam.py \
#   --bam $outPut/STARsolo_Aligned.sortedByCoord.out.bam \
#   --out_tsv $outPut/sequencing_saturation.tsv \
#   --out_png $outPut/sequencing_saturation.png \
#   --no_gene










