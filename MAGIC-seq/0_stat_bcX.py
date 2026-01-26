# 2026.1.23
# panwei
# 这个代码是用来检测MAGIC-seq数据中提取的barcode是否在白名单中, 可以只看X或者XY, 或者XYZ;
# 使用fastq.gz的文件是reformat之后的R1.fastq.gz;
# 如果有输出文件, 可以把正确的read1的名称输出出来, 输出到一个txt中;

# 脚本命令
# 只看XYZ
# python 0_stat_bcX.py --seq /data/database/MAGIC-seq-NG/20260121_second/20260121_second/00.mergeRawFq/STx170y170z7micebrain/STx170y170z7micebrain/out/STx170y170z7micebrain_reformat_R1.fastq.gz --bcx /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeA150.txt --bcy /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeB150.txt --bcz /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeC18.txt -m 3
# 只看XY
# python 0_stat_bcX.py --seq /data/database/MAGIC-seq-NG/20260121_second/20260121_second/00.mergeRawFq/STx170y170z7micebrain/Yes_haveT/T10/out/T10_reformat_R1.fastq.gz --bcx /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeA150.txt --bcy /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeB150.txt --bcz /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeC18.txt -m 2
# 只看X
# python 0_stat_bcX.py --seq /data/database/MAGIC-seq-NG/20260121_second/20260121_second/00.mergeRawFq/STx170y170z7micebrain/Yes_haveT/T10/out/T10_reformat_R1.fastq.gz --bcx /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeA150.txt --bcy /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeB150.txt --bcz /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeC18.txt -m 1
# /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeA150.txt
# /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeB150.txt
# /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeC18.txt

import argparse
import pickle
import pyfastx

# 创建参数解析器
parser = argparse.ArgumentParser(description="generate reference sequence!")
# 添加参数
parser.add_argument("--seq", required=True, help="fastq.gz")
parser.add_argument("--bcx", default=None, help="input barcodeX.txt")
parser.add_argument("--bcy", default=None, help="input barcodeY.txt")
parser.add_argument("--bcz", default=None, help="input barcodeZ.txt")
parser.add_argument("-m", "--mode", type=int, default=1, help="mode 1:barcode X, mode 2:barcode XY, mode 3:barcode XYZ")
parser.add_argument("-o", "--output", default=None, help="output.txt")

# 解析参数
args = parser.parse_args()
MODE_LEN = args.mode*8

lib_barcode_X = set()
lib_barcode_Y = set()
lib_barcode_Z = set()

# 读取barcode X 和 barcode Y, barcode Z;
if len(args.bcx) != 0:
    with open(args.bcx, 'r') as f:
        for line in f:
            # 去掉行末换行符;
            line = line.rstrip().split()
            name = line[0]
            barcode = line[1]
            lib_barcode_X.add(barcode)

if len(args.bcy) != 0:
    with open(args.bcy, 'r') as v:
        for line in v:
            line = line.rstrip().split()
            name = line[0]
            barcode = line[1]
            lib_barcode_Y.add(barcode)

if len(args.bcz) != 0:
    with open(args.bcz, 'r') as v:
        for line in v:
            line = line.rstrip().split()
            name = line[0]
            barcode = line[1]
            lib_barcode_Z.add(barcode)

Allcount = 0
count = 0
fastq = pyfastx.Fastq(args.seq)

if len(args.output) != 0:
    with open(args.output, "w", buffering=8192*1024) as f_out:
        for read in fastq:
            name = read.name
            seq = read.seq
            Allcount += 1
            if args.mode == 1:
                barcode1 = seq[:MODE_LEN]
                if barcode1 in lib_barcode_X:
                    count += 1
            elif args.mode == 2:
                barcode1 = seq[0:8]
                barcode2 = seq[8:MODE_LEN]
                if barcode1 in lib_barcode_X and barcode2 in lib_barcode_Y:
                    count += 1
            elif args.mode == 3:
                barcode1 = seq[0:8]
                barcode2 = seq[8:16]
                barcode3 = seq[16:MODE_LEN]      
                if barcode1 in lib_barcode_X and barcode2 in lib_barcode_Y and barcode3 in lib_barcode_Z:
                    count += 1
                    f_out.write(f"{name}\n")
else:
    for read in fastq:
        name = read.name
        seq = read.seq
        Allcount += 1
        if args.mode == 1:
            barcode1 = seq[:MODE_LEN]
            if barcode1 in lib_barcode_X:
                count += 1
        elif args.mode == 2:
            barcode1 = seq[0:8]
            barcode2 = seq[8:MODE_LEN]
            if barcode1 in lib_barcode_X and barcode2 in lib_barcode_Y:
                count += 1
        elif args.mode == 3:
            barcode1 = seq[0:8]
            barcode2 = seq[8:16]
            barcode3 = seq[16:MODE_LEN]      
            if barcode1 in lib_barcode_X and barcode2 in lib_barcode_Y and barcode3 in lib_barcode_Z:
                count += 1    

print("一共的reads数: ", Allcount)
print("符合", args.mode, "个barcode的reads数: ", count)
seqRatio = count/Allcount
print(f"正确barcode的比例:  {seqRatio * 100:.2f}%")

