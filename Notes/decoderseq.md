## Content
- [decoder-seq脚本解析](#decoder-seq脚本解析)




## decoder-seq脚本解析

必须准备的文件：
- barcode_coordinate.txt: 记录了barcode的顺序，比如从X1到X75和从Y1到Y75，以及对应的序列。必须包含三列，且列名称为"barcode"、"x"、"y"。"barcode"列表示所有X和Y的组合序列(16bp)，"x"列表示x轴上排序的位置，示例芯片总共75个barcode x，因此芯片x坐标从1到75。 "y"表示y轴上排序的位置，芯片y坐标从1到75。
- 全分辨率的H&E染色图。

分析流程主要包括三步：
- 1.根据barcode清单提取原始测序数据中有效的barcode。
- 2.将有效barcode的reads比对到参考基因组。
- 3.手动提供三个坐标推算每个芯片坐标的像素坐标，所有信息组织成AnnData，下游分析。

首先第1步，校正1bp错误的barcode，我的脚本只保留了1bp的错配的barcode。思路是将所有的barcode序列中的8个位置分别替换成"AGCTN"，构建真实barcode与所有候选barcode的哈希表，首先判定是否是正确的barcode，然后判定是否是候选1bp错误内的barcode。不满足条件的reads直接丢弃。最后使用pigz压缩为gz文件。
```shell
FQ1=DCseqformalMB_S1_L001_R1_001.fastq.gz
FQ2=DCseqformalMB_S1_L001_R2_001.fastq.gz
barcode=7.5mm_barcode_coordinate.txt
outPut=./result
Numthread=16
# 脚本在Decoder-seq文件夹中
python script1_getBarcodeSR_decoderseq.py \
        -i $FQ1 -I $FQ2 \
        -x $barcode \
        -o $outPut
pigz -p $Numthread $outPut/DCseqformalMB_S1_L001_R1_001_trim.fastq
pigz -p $Numthread $outPut/DCseqformalMB_S1_L001_R2_001_trim.fastq
```
其次第二步，使用[STARsolo](https://github.com/weiwei4396/Spatial-Transcriptomics/blob/main/Notes/annData.md) 将过滤后的reads比对到参考基因组。
```shell
sampath=./result
sample=DCseqformalMB_S1_L001
fastq1=${sampath}/${sample}_R1_001_trim.fastq.gz
fastq2=${sampath}/${sample}_R2_001_trim.fastq.gz
# 这个barcode只需要一列，就是所有barcode的组合，假设barcode x有75个, barcode y有75个，那么这个文件就是5625行1列的文本文件。
whitelist=./result/second_column.txt
# 注释文件
MAP=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/star2710b
ANN=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/genes/genes.gtf
t_num=32
ulimit -n 100000
mkdir ${sampath}/STARsolo

/data/workdir/panw/software/STAR-2.7.11b/bin/Linux_x86_64/STAR --genomeDir ${MAP} \
  --outFileNamePrefix ${sampath}/STARsolo/${sample}_ \
  --readFilesCommand zcat \
  --readFilesIn ${fastq2} ${fastq1} \
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
最后第3步，从STARsolo读取每个barcode/spot的基因表达矩阵；再从barcode_coordinate.txt读取每个barcode对应的芯片网格坐标(x, y)；根据用户在H&E图像上手动提供的三个角点(初始芯片位置的像素坐标)像素坐标，线性推算出整个75×75芯片网格中每个网格坐标对应的H&E图像像素坐标，这样得到了芯片坐标和像素坐标的哈希表；然后把每个barcode的表达矩阵、芯片坐标、像素坐标、组织区域信息和H&E图像一起组织成 AnnData，用于后续空间可视化、QC 和分析。
```shell
python script3_decoderseq_downstream.py
```
需要改动的参数如下，STARsolo的结果必须将filtered中的matrix.mtx改为matrix.mtx.gz
```shell
# SAMPLE = 'sample1'
# STARSOLO_DIR = '/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/bingqi_Decoderseq/DCseqformalMB_S1_L001_Solo.out'
# BARCODE_POS  = '/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/bingqi_Decoderseq/7.5mm_barcode_coordinate.txt'
# HE_IMAGE     = '/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/bingqi_Decoderseq/MyflipTest.tif'
# outPut       = '/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/bingqi_Decoderseq/pictures'
# os.makedirs(f'{outPut}', exist_ok=True)
# 这里给的坐标和barcode不是很清楚, 通过最后结果判断是翻转了;
# 默认不做操作的话用identity, 这个示例数据需要翻转一下;
# coord_transform_mode = 'flip_y' # identity/swap_xy/rotate_ccw_90/rotate_cw_90/...
# ----- Step 1: get the three corner pixel positions ---------------
# 获取左上, 右上, 左下的三个像素的位置;
# corner_X1Y1  = (480, 464)    # top-left
# corner_X75Y1 = (16144, 496)  # top-right
# corner_X1Y75 = (528, 16128)  # bottom-left
```







