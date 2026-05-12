## Content
- [MAGIC-seq脚本解析](#MAGIC-seq脚本解析)

## **MAGIC-seq脚本解析**
### 我自己写的步骤
必须准备的文件：两个/三个barcode白名单文件
- [Spatial_barcodeA150.txt](https://github.com/weiwei4396/Spatial-Transcriptomics/blob/main/picture/magicseq_whitelist.png)：沿用的MAGIC-seq文章数据给定的白名单结构，第一列表示芯片上的位置，第二列表示这个barcode的序列。
- Spatial_barcodeB150.txt：跟barcode x情况类似，且数量相同。
- Spatial_barcodeC18.txt(看是否有Z轴)：情况也相似，只不过数量很少。
- 全分辨率的H&E染色图。

分析流程主要包括三步：
- 1.根据barcode清单提取原始测序数据中有效的barcode。[script1_getBarcodeSR_magicseq.py](https://github.com/weiwei4396/Spatial-Transcriptomics/blob/main/MAGIC-seq/script1_getBarcodeSR_magicseq.py)
- 2.将有效barcode的reads比对到参考基因组。
- 3.手动提供三个坐标推算每个芯片坐标的像素坐标，所有信息组织成AnnData，下游分析。

首先第1步，校正1bp错误的barcode，我的脚本只保留了1bp的错配的barcode。思路是将所有的barcode序列中的8个位置分别替换成"AGCTN"，构建真实barcode与所有候选barcode的哈希表，首先判定是否是正确的barcode，然后判定是否是候选1bp错误内的barcode。不满足条件的reads直接丢弃。最后使用pigz压缩为gz文件。

```shell
FQ1=/data/database/MAGIC-seq-NG/E17-1/CRR1158992_R1.fastq.gz
FQ2=/data/database/MAGIC-seq-NG/E17-1/CRR1158992_R2.fastq.gz
barcodeX=/data/database/MAGIC-seq-NG/E17-1/Barcode-M9-150-E17.5/Spatial_barcodeA150.txt
barcodeY=/data/database/MAGIC-seq-NG/E17-1/Barcode-M9-150-E17.5/Spatial_barcodeB150.txt
barcodeZ=/data/database/MAGIC-seq-NG/E17-1/Barcode-M9-150-E17.5/Spatial_barcodeC18.txt
resultLinker=/data/database/MAGIC-seq-NG/E17-1/result
getPY=script1_getBarcodeSR_magicseq.py

python $getPY \
        -i $FQ1 -I $FQ2 \
        --bcx $barcodeX --bcy $barcodeY --bcz $barcodeZ \
        -m 3 -o $resultLinker
cd result
pigz -p 16 CRR1158992_R1_trim.fastq
pigz -p 16 CRR1158992_R2_trim.fastq
```
其次第二步，使用[STARsolo](https://github.com/weiwei4396/Spatial-Transcriptomics/blob/main/Notes/annData.md) 将过滤后的reads比对到参考基因组。
```shell
# 1.准备数据
sampath=/data/database/MAGIC-seq-NG/E17-1/resultYihuan
sample=CRR1158992
fastq1=${sampath}/${sample}_R1_trim.fastq.gz
fastq2=${sampath}/${sample}_R2_trim.fastq.gz
# 给定的白名单
whitelist=/data/database/MAGIC-seq-NG/E17-1/Barcode-M9-150-E17.5/whitelist.txt
#ID=/data/database/MAGIC-seq-NG/E17-1/Barcode-M9-150-E17.5/M9_ST_ids_barcode_chip1_C18.txt
# 参考基因组
MAP=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/star2710b
ANN=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/genes/genes.gtf
# 线程
t_num=16
ulimit -n 100000
mkdir ${sampath}/STARsolo
# 2.比对
/data/workdir/panw/software/STAR-2.7.11b/bin/Linux_x86_64/STAR --genomeDir ${MAP} \
  --outFileNamePrefix ${sampath}/STARsolo/${sample}_ \
  --readFilesCommand pigz -p 8 -dc \
  --readFilesIn ${fastq2} ${fastq1} \
  --outSAMattributes NH HI nM AS CR UR CY UY CB UB GX GN sS sQ sM sF \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 121539607552 \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist ${whitelist} \
  --soloCBstart 1 \
  --soloCBlen 24 \
  --soloUMIstart 25 \
  --soloUMIlen 12 \
  --soloFeatures Gene GeneFull SJ Velocyto \
  --soloMultiMappers EM \
  --soloUMIdedup 1MM_All \
  --soloCellFilter EmptyDrops_CR \
  --soloCellReadStats Standard \
  --clipAdapterType CellRanger4 \
  --outReadsUnmapped Fastx \
  --runThreadN ${t_num}
```
- 上面的代码是针对的三个barcode的情况，两个barcode的情况同理，可以直接看脚本：[BarcodeXY](https://github.com/weiwei4396/Spatial-Transcriptomics/blob/main/MAGIC-seq/BarcodeXY.sh) | [BarcodeXYZ](https://github.com/weiwei4396/Spatial-Transcriptomics/blob/main/MAGIC-seq/BarcodeXYZ.sh)
- 额外步骤, 提取前5000条序列blastn查看是否有污染
```python
# 脚本/data/database/MAGIC-seq-NG/Spatial-Transcriptomics/MAGIC-seq/1_blast_species.py
# 给定read2.fastq.gz可以自动提取前5000条序列, 也可以直接给定fasta或者fastq
# 最后需要给定nt库
python 1_blast_species.py -q RNA20X125Y1_raw_R2.fastq.gz -d /data/workdir/zhangj/database/nt_/core/core_nt
```
最后第3步，从STARsolo读取每个barcode/spot的基因表达矩阵；再从barcode_coordinate.txt读取每个barcode对应的芯片网格坐标(x, y)；根据用户在H&E图像上手动提供的三个角点(初始芯片位置的像素坐标)像素坐标，线性推算出整个75×75芯片网格中每个网格坐标对应的H&E图像像素坐标，这样得到了芯片坐标和像素坐标的哈希表；然后把每个barcode的表达矩阵、芯片坐标、像素坐标、组织区域信息和H&E图像一起组织成 AnnData，用于后续空间可视化、QC 和分析。

___

### MAGIC-seq文章的pipline

BarcodeX + BarcodeY
<details>
<summary> X+Y barcode流程 </summary>

- 数据预处理, 提取Read1的barcode和UMI, 得到只由Barcode和UMI组成序列; 根据实验设计是否有barcodeZ, 提取不同的索引; 三宫格或九宫格: BarcodeX[1-8], BarcodeY+UMI[27-46];
- 整个切片: BarcodeX[1-8], BarcodeY[27-34], BarcodeZ+UMI[53:72];
- STARsolo所需的read1格式为barcode+UMI; 因此设置为**BarcodeX/BarcodeY/UMI**或**BarcodeX/BarcodeY/BarcodeZ/UMI**;

1.数据准备
```shell
# 设置输入文件; sampath为当前测序的R1.fastq.gz和R2.fastq.gz所在的文件夹;
sampath=/data/database/MAGIC-seq-NG/Olfb
sample=OlfBulb
fastq1=${sampath}/${sample}_R1.fastq.gz
fastq2=${sampath}/${sample}_R2.fastq.gz

# whitelist为所有的barcode的组合, 如X和Y共4900种barcode, 每个barcode占一行;
whitelist=/data/database/MAGIC-seq-NG/Olfb/Mouse_Adult_Organ_T9_70_50um/whitelist
ID=/data/database/MAGIC-seq-NG/Olfb/Mouse_Adult_Organ_T9_70_50um/T9-ids-barcode.txt

# 参考基因组, 小鼠/人;
# 创建STAR的参考基因组index;
# /data/workdir/panw/software/STAR-2.7.11b/bin/Linux_x86_64/STAR --runThreadN 16 --runMode genomeGenerate --genomeDir star2710b --genomeFastaFiles /data/workdir/panw/reference/human/refdata-gex-GRCh38-2024-A/fasta/genome.fa --sjdbGTFfile /data/workdir/panw/reference/human/refdata-gex-GRCh38-2024-A/genes/genes.gtf --sjdbOverhang 100 --genomeSAindexNbases 14 --genomeChrBinNbits 18 --genomeSAsparseD 3
# MAP=/data/workdir/panw/reference/human/refdata-gex-GRCh38-2024-A/star2710b
# ANN=/data/workdir/panw/reference/human/refdata-gex-GRCh38-2024-A/genes/genes.gtf
MAP=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/star2710b
ANN=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/genes/genes.gtf

# 线程;
t_num=16

# 建立输出的文件和文件夹;
log_file=${sampath}/${sample}_st_log.txt
st_path=${sampath}/${sample}
touch log_file
mkdir ${st_path}
mkdir ${st_path}/split
mkdir ${st_path}/out
```

2.提取barcode和umi
```shell
# 注意: seqkit需要版本 2.0.0, 至少已知 2.10.0 中 concat功能不能支持现在的分析;
# 2.seqkit分割文件处理; 默认分成10份; 然后分别提取 BarcodeX和BarcodeY-UMI;

FileNum=10
seqkit split2 -1 ${sampath}/${sample}_R1.fastq.gz -2 ${sampath}/${sample}_R2.fastq.gz -p $FileNum -O ${st_path}/split -f

files=(001 002 003 004 005 006 007 008 009 010)

for file in "${files[@]}"
do
    seqkit subseq -r 1:8 ${st_path}/split/${sample}_R1.part_${file}.fastq.gz -o ${st_path}/out/part_${file}_test1-8.fastq.gz
    seqkit subseq -r 27:46 ${st_path}/split/${sample}_R1.part_${file}.fastq.gz -o ${st_path}/out/part_${file}_test27-46.fastq.gz
    seqkit concat ${st_path}/out/part_${file}_test1-8.fastq.gz ${st_path}/out/part_${file}_test27-46.fastq.gz -o ${st_path}/out/${sample}_${file}_reformat_R1.fastq.gz
    seqkit replace -p " (.*)$" -r "" ${st_path}/split/${sample}_R2.part_${file}.fastq.gz -o ${st_path}/out/${sample}_${file}_reformat_R2.fastq.gz
    rm -rf ${st_path}/out/part_${file}_test1-8.fastq.gz
    rm -rf ${st_path}/out/part_${file}_test27-46.fastq.gz
done


cat ${st_path}/out/*_reformat_R1.fastq.gz > ${st_path}/out/${sample}_reformat_R1.fastq.gz
rm -rf ${st_path}/out/${sample}_0*_reformat_R1.fastq.gz
cat ${st_path}/out/*_reformat_R2.fastq.gz > ${st_path}/out/${sample}_reformat_R2.fastq.gz
rm -rf ${st_path}/out/${sample}_0*_reformat_R2.fastq.gz
rm -rf ${st_path}/split
```

3.质控
```shell
# TA建库需要运行这一步去接头;
# AAGCAGTGGTATCAACGCAGAGTGAATGGG
cutadapt -g AAGCAGTGGTATCAACGCAGAGTGAATGGG -e 0.01 -j ${t_num} -o ${st_path}/${sample}_reformat_cutadapt_R2.fastq.gz ${st_path}/out/${sample}_reformat_R2.fastq.gz

fastp -i ${st_path}/out/${sample}_reformat_R1.fastq.gz -I ${st_path}/${sample}_reformat_cutadapt_R2.fastq.gz \
       -o ${st_path}/${sample}_cleaned_R1.fastq.gz -O ${st_path}/${sample}_cleaned_R2.fastq.gz \
       -l 28 -x -g -w ${t_num} --detect_adapter_for_pe -j ${st_path}/${sample}.barcode.fastp.json -h ${st_path}/${sample}.barcode.fastp.html
```


4.比对; 将数据比对到参考基因组, 需要注意的是提供的reads先输入reads2再输入reads1;
```shell
# 使用STARsolo将数据比对;
/data/workdir/panw/software/STAR-2.7.11b/bin/Linux_x86_64/STAR --genomeDir ${MAP} \
  --outFileNamePrefix ${st_path}/STARsolo/${sample}_ \
  --readFilesCommand zcat \
  --readFilesIn ${st_path}/${sample}_cleaned_R2.fastq.gz ${st_path}/${sample}_cleaned_R1.fastq.gz \
  --outSAMattributes NH HI nM AS CR UR CY UY CB UB GX GN sS sQ sM sF \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 121539607552 \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist ${whitelist} \
  --soloCBstart 1 \
  --soloCBlen 16 \
  --soloUMIstart 17 \
  --soloUMIlen 12 \
  --soloFeatures Gene GeneFull SJ Velocyto \
  --soloMultiMappers EM \
  --soloUMIdedup 1MM_All \
  --soloCellFilter EmptyDrops_CR \
  --soloCellReadStats Standard \
  --clipAdapterType CellRanger4 \
  --outReadsUnmapped Fastx \
  --runThreadN ${t_num}
```

5.将STARsolo生成的基因表达转化为Anndata;

</details>


BarcodeX + BarcodeY + BarcodeZ
<details>
<summary> X+Y+Z barcode流程 </summary>

- 数据预处理, 提取Read1的barcode和UMI, 得到只由Barcode和UMI组成序列; 根据实验设计是否有barcodeZ, 提取不同的索引; 三宫格或九宫格: BarcodeX[1-8], BarcodeY+UMI[27-46];
- 整个切片: BarcodeX[1-8], BarcodeY[27-34], BarcodeZ+UMI[53:72];
- STARsolo所需的read1格式为barcode+UMI; 因此设置为**BarcodeX/BarcodeY/UMI**或**BarcodeX/BarcodeY/BarcodeZ/UMI**;

1.数据准备
```shell
# 设置输入文件; sampath为当前测序的R1.fastq.gz和R2.fastq.gz所在的文件夹;
sampath=/data/database/MAGIC-seq-NG/20260121_second/20260121_second/00.mergeRawFq/STx170y170z7micebrain
sample=STx170y170z7micebrain
fastq1=${sampath}/${sample}_R1.fastq.gz
fastq2=${sampath}/${sample}_R2.fastq.gz

# whitelist为所有的barcode的组合, 如X和Y共4900种barcode, 每个barcode占一行;
whitelist=/data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/whitelist.txt
ID=/data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/M9_ST_ids_barcode_chip1_C18.txt

# 参考基因组, 小鼠/人;
# 创建STAR的参考基因组index;
# /data/workdir/panw/software/STAR-2.7.11b/bin/Linux_x86_64/STAR --runThreadN 16 --runMode genomeGenerate --genomeDir star2710b --genomeFastaFiles /data/workdir/panw/reference/human/refdata-gex-GRCh38-2024-A/fasta/genome.fa --sjdbGTFfile /data/workdir/panw/reference/human/refdata-gex-GRCh38-2024-A/genes/genes.gtf --sjdbOverhang 100 --genomeSAindexNbases 14 --genomeChrBinNbits 18 --genomeSAsparseD 3
# MAP=/data/workdir/panw/reference/human/refdata-gex-GRCh38-2024-A/star2710b
# ANN=/data/workdir/panw/reference/human/refdata-gex-GRCh38-2024-A/genes/genes.gtf
MAP=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/star2710b
ANN=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/genes/genes.gtf

# 线程;
t_num=16

# 建立输出的文件和文件夹;
log_file=${sampath}/${sample}_st_log.txt
st_path=${sampath}/${sample}
touch log_file
mkdir ${st_path}
mkdir ${st_path}/split
mkdir ${st_path}/out
```

</details>

### 之前自写的对MAGIC-seq数据排查的脚本

<details>
<summary> MAGIC-seq数据中对建库问题的排查 </summary>

```shell
conda activate py310
```

0.通过每个位置碱基含量百分比, 初步判断read1和read2的结果;
```shell
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
```

1.按照确定位置找barcode, 查看正确的barcode比例; 通过改变-m参数可以查看m个barcode是否正确;

```python
python 0_stat_bcX.py --seq /data/database/MAGIC-seq-NG/20260121_second/20260121_second/00.mergeRawFq/STx170y170z7micebrain/Yes_haveT/T10/out/T10_reformat_R1.fastq.gz --bcx /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeA150.txt --bcy /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeB150.txt --bcz /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeC18.txt -m 1
```

2.查看read2中TSO含量的比例;
```python
python 2_TSO_ratio.py --tso AAGCAGTGGTATCAACGCAG --fq STx170y170z7micebrain_R2.fastq.gz
```

3.抽取部分reads使用blastn比对查看一些物种信息;
```shell
fa=/data/database/MAGIC-seq-NG/20260121_second/20260121_second/00.mergeRawFq/STx170y170z7micebrain/STAll_R2.fasta
blastn -query ${fa} -db /data/workdir/zhangj/database/nt_/core/core_nt \
        -out ./blastn_results.tsv -max_target_seqs 5 -task blastn \
        -outfmt "6 std qlen slen btop ssciname staxid stitle" \
        -num_threads 30
```
使用python脚本解析blastn_results.tsv中包含的物种信息;

4.检查read1中barcode的自连情况, 顺便也可以看看barcode的连接情况;
```python
python barcodeXXX.py --fq STx170y170z7micebrain_R1.fastq.gz --bcx /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeA150.txt --bcy /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeB150.txt --bcz /data/database/MAGIC-seq-NG/P0-1/Barcode-M9-150-P04/Spatial_barcodeC18.txt -m 3
```

</details>






