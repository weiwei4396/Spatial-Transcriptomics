# 2026.1.26
# panwei
# 这个代码是用来查看TSO在read2中占比;
# 脚本命令
# python 2_TSO_ratio.py -t AAGCAGTGGTATCAACGCAGAGTGAATGGG -q oneBCST293TRNA_2.fq.gz

import argparse
import pyfastx

# 创建参数解析器
parser = argparse.ArgumentParser(description="generate reference sequence!")
# 添加参数
parser.add_argument("-t", "--tso", required=True, help="which TSO")
parser.add_argument("-q", "--fq", required=True, help="read2.fq.gz")

# 解析参数
args = parser.parse_args()

Allcount = 0
TSO_arise = 0
# TSO_count = defaultdict(int) 
TSO_count = dict()

fastq = pyfastx.Fastq(args.fq)
for read in fastq:
    Allcount += 1
    name = read.name
    seq = read.seq
    if Allcount % 1000000 == 0:
        print(f"已经处理了{Allcount/1000000}M reads...")
    if args.tso in seq:
        TSO_arise += 1
        count = seq.count(args.tso)
        TSO_count[count] = TSO_count.get(count, 0) + 1

print("要看的TSO序列是:", args.tso, "[AAGCAGTGGTATCAACGCAGAGTGAATGGG]")
print(f"总共的reads数: {Allcount}")
arise_ratio = TSO_arise/Allcount
print(f"TSO出现在 {TSO_arise} 条reads中, 占比{arise_ratio:.2%}")

for eachcondition, eachcount in TSO_count.items():
    arise_ratio = eachcount/Allcount
    print(f"TSO在{eachcount}条reads中出现了{eachcondition}次, 占比{arise_ratio:.2%};")





















