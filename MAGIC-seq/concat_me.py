# 2025.11.24
# panwei
# 这个代码是用来替代MAGIC-seq流程中合并seqkit concat两个文件的功能;
# 脚本命令


import argparse
import re
import gzip

# 创建参数解析器
parser = argparse.ArgumentParser(description="seqkit concat !")
# 添加参数
parser.add_argument("-a", "--partA", required=True, help="concat 1-8 barcode")
parser.add_argument("-b", "--partB", required=True, help="concat 27-46 barcode and umi")
parser.add_argument("-o", "--output", required=True, help="output concat files")

# 解析参数
args = parser.parse_args()

last_twoA = args.partA[-2:]
last_twoB = args.partB[-2:]

row_count = 0

if last_twoA == "gz" and last_twoB == "gz":
    with gzip.open(args.partA, "rt") as fa, gzip.open(args.partB, "rt") as fb, gzip.open(args.output, "wt") as outfile:
        for lineA, lineB in zip(fa, fb):
            row_count = row_count + 1
            if row_count%4 == 1:
                outfile.write(lineA)
            elif row_count%4 == 2:
                lineA = lineA.strip()
                lineB = lineB.strip()
                outfile.write(lineA+lineB)
                outfile.write("\n")
            elif row_count%4 == 3:
                outfile.write(lineA)
            else:
                lineA = lineA.strip()
                lineB = lineB.strip()
                outfile.write(lineA+lineB)
                outfile.write("\n")           

else :
    with open(args.partA, "r") as fa, open(args.partB, "r") as fb, open(args.output, "w") as outfile:
        for lineA, lineB in zip(fa, fb):
            row_count = row_count + 1
            if row_count%4 == 1:
                outfile.write(lineA)
            elif row_count%4 == 2:
                lineA = lineA.strip()
                lineB = lineB.strip()
                outfile.write(lineA+lineB)
                outfile.write("\n")
            elif row_count%4 == 3:
                outfile.write(lineA)
            else:
                lineA = lineA.strip()
                lineB = lineB.strip()
                outfile.write(lineA+lineB)
                outfile.write("\n")                


