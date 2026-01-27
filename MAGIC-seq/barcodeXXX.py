# 2026.01.27
# panwei
# 这个脚本是想检查一下barcode的自连情况;
# python barcodeXXX.py --fq STx170y170z7micebrain_R1.fastq.gz --bcx /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeA150.txt --bcy /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeB150.txt --bcz /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeC18.txt -m 3

import argparse
import pyfastx

# 创建参数解析器
parser = argparse.ArgumentParser(description="generate reference sequence!")
# 添加参数
parser.add_argument("-q", "--fq", required=True, help="input R1 fastq.gz")
parser.add_argument("-x", "--bcx", default=None, help="input barcodeX whitelist")
parser.add_argument("-y", "--bcy", default=None, help="input barcodeY whitelist")
parser.add_argument("-z", "--bcz", default=None, help="input barcodeZ whitelist")
parser.add_argument("-m", "--mode", default=3, help="查看哪个barcode的自连, 1 X 2 Y 3 Z")

# 解析参数
args = parser.parse_args()

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


fastq = pyfastx.Fastq(args.fq)

Allcount = 0
Xcount = 0
XXcount = 0
XYcount = 0
XZcount = 0
XLinker = 0
XLYcount = 0
XLXcount = 0
XLZcount = 0
XLLcount = 0
XLYLcount = 0
XLYLZcount = 0

for read in fastq:
    name = read.name
    seq = read.seq
    Allcount += 1
    if Allcount % 1000000 == 0:
        print(f"have processed {Allcount/1000000}M reads...")
    
    barcodeX = seq[:8]
    if barcodeX in lib_barcode_X:
        Xcount += 1
        barcodeXX = seq[8:16]
        if barcodeXX in lib_barcode_X:
            XXcount += 1
        elif barcodeXX in lib_barcode_Y:
            XYcount += 1
        elif barcodeXX in lib_barcode_Z:
            XZcount += 1
        elif barcodeXX == "CAGTCATG":
            XLinker += 1
            barcodeYY = seq[26:34]
            if barcodeYY in lib_barcode_Y:
                XLYcount += 1
                if args.mode == 3:
                    LinkerYZ = seq[34:42]
                    if LinkerYZ == "TGATGCGA":
                        XLYLcount += 1
                        barcodeZZ = seq[52:60]
                        if barcodeZZ in lib_barcode_Z:
                            XLYLZcount += 1
            elif barcodeYY in lib_barcode_X:
                XLXcount += 1
            elif barcodeYY in lib_barcode_Z:
                XLZcount += 1
            elif barcodeYY == "TGATGCGA":
                XLLcount += 1


print(f"reads总数:{Allcount}")
print(f"reads中barcode X连接正确的总数:{Xcount}")
print(f"reads中barcode X后面紧跟另一个X的总数:{XXcount}")
print(f"reads中barcode X后面紧跟Y的总数:{XYcount}")
print(f"reads中barcode X后面紧跟Z的总数:{XZcount}")
print(f"reads中barcode X后面紧跟Linker的总数:{XLinker}")
print(f"reads中barcode X后面紧跟Linker又连对了barcode Y总数:{XLYcount}")
print(f"reads中barcode X后面紧跟Linker又连了barcode X总数:{XLXcount}")
print(f"reads中barcode X后面紧跟Linker又连了barcode Z总数:{XLZcount}")
print(f"reads中barcode X后面紧跟LinkerXY又连了LinkerYZ总数:{XLLcount}")
print(f"reads中barcode X后面紧跟LinkerXY又连对了barcode Y又连上了LinkerYZ总数:{XLYLcount}")
print(f"reads中barcode X后面紧跟LinkerXY又连对了barcode Y又连上了LinkerYZ又连上了barcode Z总数:{XLYLZcount}")

