## decoder-seq脚本

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
其次第二步，使用STARsolo(https://github.com/weiwei4396/Spatial-Transcriptomics/blob/main/Notes/annData.md) 将过滤后的reads比对到参考基因组。

最后第3步，从 STARsolo 读取每个 barcode/spot 的基因表达矩阵；再从 barcode_coordinate.txt 读取每个 barcode 对应的芯片网格坐标 (x, y)；根据用户在 H&E 图像上手动提供的三个角点(初始芯片位置的像素坐标)像素坐标，线性推算出整个 75×75 芯片网格中每个网格坐标对应的 H&E 图像像素坐标，这样得到了芯片坐标和像素坐标的哈希表；然后把每个 barcode 的表达矩阵、芯片坐标、像素坐标、组织区域信息和 H&E 图像一起组织成 AnnData，用于后续空间可视化、QC 和分析。







