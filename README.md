# Spatial-Transcriptomics
Learning some analysis code.

## 0.**Introduction**
### Sequencing-based spatial transcriptomics
### Imaging-based spatial transcriptomics

基于成像的平台（也称为基于单分子检测的平台）通过原位杂交（ISH）或原位测序（ISS）技术，对预设的靶向基因面板（通常包含数百至数千个基因）进行逐轮检测，从而在组织中原位识别单个RNA分子的空间位置。由于每个转录本都被单独识别，原始数据天然具有亚细胞级别的空间分辨率。

在预处理阶段，需通过图像分割来识别单个细胞或细胞核的边界，并将检测到的 RNA分子分配给相应的细胞或细胞核。然而，细胞分割极具挑战性，尤其在组织切片中：由于细胞在三维空间中存在重叠，在二维切片的某一 (x, y) 位置可能包含多个细胞的信号，导致边界难以准确界定。

完成分割后，可将RNA分子按细胞汇总为细胞级别的基因表达计数，也可直接在单分子水平开展分析。若采用细胞级别数据，许多分析方法可复用自基于“spot”的空间转录组学或单细胞 RNA 测序领域已有的工具。



## 1.**SpaceRanger**
[官方网址](https://www.10xgenomics.com/support/software/space-ranger/latest)



## 2.**AnnData data structure**

[AnnData结构](https://raw.githubusercontent.com/weiwei4396/Spatial-Transcriptomics/main/picture/anndata.jpg)
scverse 是一个专注于生命科学基础工具的组织和生态系统，最初聚焦于单细胞数据分析。它的优势在于出色的扩展性、灵活性以及与现有Python数据科学和机器学习工具的强大互操作性。

在scverse生态系统中, AnnData是用来将数据矩阵与这些注释关联起来的核心工具。为了提高效率, AnnData支持稀疏矩阵 (sparse matrices) 和部分读取(partial reading), 这样可以更快地处理大规模数据。AnnData在功能上与R生态系统中的数据结构 (比如Bioconductor的SummarizedExperiment或Seurat对象)相似, 但R包通常使用转置后的特征矩阵 (基因 x 细胞)。

在AnnData的核心中, 存储了一个稀疏或密集矩阵 (在scRNA-seq中就是计数矩阵), 称为X, 这个矩阵的维度是 obs_names x var_names (细胞 x 基因), 其中obs(观测值)对应细胞条形码, var(变量)对应基因标识符。 矩阵X被两个Pandas数据框(DataFrame)包围。其中obs保存细胞的注释信息, var保存基因的注释信息。 

AnnData还可以储存很多额外信息. 比如, 其他关于观测值和变量的多维数据(如UMAP) 储存在obsm和varm中, 图结构(比如细胞之间的关系或基因之间的关系)存储在obsp和varp中, 任何不适合其他槽位的非结构化数据都可以存储在uns中, 还可以通过layers存储矩阵X的额外值。例如, 可以在名为counts的层中存储未经标准化的原始计数数据, 而在默认的层中存储标准化后的数据。

[AnnData](https://mp.weixin.qq.com/s/0OFRSB3BZcNltHkp_1VG1Q)

<details>
<summary> </summary>

- obs: 'in_tissue', 'array_row', 'array_col', 'gender', 'age', 'tissue', 'Strain', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'n_counts', 'n_genes'
- var: 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
- uns: 'spatial', 'log1p', 'hvg', 'pca', 'neighbors', 'umap', 'tsne', 'tissue_colors'
- obsm: 'spatial', 'X_pca', 'X_umap', 'X_tsne'
- varm: 'PCs'
- obsp: 'distances', 'connectivities'

</details>

## 3.**SpaceRanger outs folder**

[outs输出文件结构](https://github.com/weiwei4396/Spatial-Transcriptomics/blob/main/picture/SpaceRanger_outs.jpg)

使用 [SpaceRanger](https://www.cnblogs.com/huanping/p/16839765.html)

<details>
<summary> </summary>
  
- filtered_feature_bc_matrix [folder] 跟单细胞一样的三个.gz文件，矩阵、基因名称和barcode名称。
- filtered_feature_bc_matrix.h5 本质上跟上面的文件夹存储的信息是一致的，都是空间转录组表达矩阵的信息。还包含了一些metadata信息，包括SpaceRanger版本信息，测序实验试剂信息等，直接读取这个h5文件包含的信息更多。
- spatial [folder] 这个文件夹中包含了几个文件，主要包含空转切片的图片信息。
  + detected_tissue_image.jpg：红色或者蓝色表示被测序捕获到的区域，灰色区域表示切片没有被测序捕获到的区域;
  + tissue_hires_image.png: 这两个图片都是对原始图片的下采样，方便对展示数据，不然图片太大了;
  + tissue_lowres_image.png: 更低分辨率的下采样;
  + scalefactors_json.json: 图片之间的数值转换关系;
  + tissue_positions_list.csv: 六列信息，barocode，spot是否在tissue里，spot所在的行，spot所在的列，细胞所在图像上的位置行和列(后两个);
  + aligned_tissue_image.jpg

</details>


## 4.**MAGIC-seq & Decoder-seq Analysis Pipline**

### 4.1 自写的对于MAGIC-seq和Decoder-seq的pipline

- 第一步，通过**getBarcode_SR_BroCOLI.py**获取barcode连接的信息, 提取出带有可信barcode的序列。
- 脚本中有对barcode的校正，只针对了错配错误的校正，没有针对插入和删除。通过枚举错误，在集合中哈希查找，速度相对更快。
```
# 添加write_discarded参数可以输出那些barcode不准确的reads, 默认不加这个参数。
barcodeX=/data/database/MAGIC-seq-NG/Spatial-Transcriptomics/bingqi70_X.txt
barcodeY=/data/database/MAGIC-seq-NG/Spatial-Transcriptomics/bingqi70_Y.txt
python getBarcode_SR_BroCOLI.py -i R1.fastq.gz -I R2.fastq.gz -x $barcodeX -y $barcodeY -m 2 -o $resultPath --write_discarded

# 加上z之后用法
barcodeX=/data/database/MAGIC-seq-NG/E17-1/Barcode-M9-150-E17.5/Spatial_barcodeA150.txt
barcodeY=/data/database/MAGIC-seq-NG/E17-1/Barcode-M9-150-E17.5/Spatial_barcodeB150.txt
barcodeZ=/data/database/MAGIC-seq-NG/E17-1/Barcode-M9-150-E17.5/Spatial_barcodeC18.txt
python getBarcode_SR_BroCOLI.py -i R1.fastq.gz -I R2.fastq.gz -x $barcodeX -y $barcodeY -z $barcodeZ -m 3 -o $resultPath --write_discarded
```

- 第二步, 使用STARsolo比对, 这个步骤跟magic-seq基本相同。
```shell
# 1.准备数据和文件;
sampath=/data/database/MAGIC-seq-NG/20260417_geneadd/gene0417/00.mergeRawFq/RNA20X125Y1/resultYihuan
sample=RNA20X125Y1_raw
fastq1=${sampath}/${sample}_R1_trim.fq.gz
fastq2=${sampath}/${sample}_R2_trim.fq.gz
# 给定的白名单;
whitelist=/data/database/MAGIC-seq-NG/Olfb/Mouse_Adult_Organ_T9_70_50um/whitelist
ID=/data/database/MAGIC-seq-NG/Olfb/Mouse_Adult_Organ_T9_70_50um/T9-ids-barcode.txt
MAP=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/star2710b
ANN=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/genes/genes.gtf
# 线程;
t_num=16
ulimit -n 100000
mkdir ${sampath}/STARsolo
# 2.比对
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
- 额外步骤, 提取前5000条序列blastn查看是否有污染
```python
# 脚本/data/database/MAGIC-seq-NG/Spatial-Transcriptomics/MAGIC-seq/1_blast_species.py
# 给定read2.fastq.gz可以自动提取前5000条序列, 也可以直接给定fasta或者fastq
# 最后需要给定nt库
python 1_blast_species.py -q RNA20X125Y1_raw_R2.fastq.gz -d /data/workdir/zhangj/database/nt_/core/core_nt
```

### 4.2 MAGIC-seq的pipline

BarcodeX+BarcodeY, 两种Barcode
<details>
<summary> </summary>

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


BarcodeX+BarcodeY+BarcodeZ, 三种Barcode
<details>
<summary> </summary>

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

### 4.3 之前自写的对MAGIC-seq数据排查的脚本

MAGIC-seq数据中对建库问题的排查
<details>
<summary> </summary>


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

### 4.4 decoder-seq的脚本
[Decoder-seq分析流程](https://github.com/songjiajia2018/Decoder-seq/tree/main)

测试过这个分析流程，能根据编辑距离很好的校正barcode，缺点是相对较慢，与我自己写的脚本的区别在于它能校正1bp的插入和删除的barcode错误，我放弃了这部分，只保留了1bp的错配的barcode。而MAGIC-seq更粗暴，只是按位置提取出来，只保留了正确的barcode。


## 参考

[空转的分析模式](https://mp.weixin.qq.com/s/lnP752JylTJouBwkUuGT-Q)


## **Question & Answer**
<details>
<summary> </summary>

### 1. Segment Anything Model (SAM)



### 2. image registration (图像配准)

- 图像配准就是把两张不同来源, 不同大小, 不同角度的图像对齐, 使它们的同一个位置在像素空间上对应起来;
- H&E图是经过扫描仪拍摄得到的, 随着组织样本放置, 扫描角度、形变, 不同批次 H&E 图像中的网格, 不一定是正交的、等距的;
- 在H&E图和对应的spatial spot图上的同一个位置选取参考点, 两个对应点通常用于仿射变换, 线性缩放, 旋转, 平移;
- 像素坐标: 每个像素有自己的位置, 左上角永远是 (0, 0), 向右为x+, 向下为y+;


### 3. Gene+ 数据
```shell
wget -c --no-check-certificate "http://geneplus001.oss-cn-beijing.aliyuncs.com/kefu%2F20260105_ZhongGuoKeXueYuanBeiJingJiYinZuYanJiuSuoGuoJiaShengWuXinXiZhongXin-wubingqi-1_1.tar?Expires=1769112506&OSSAccessKeyId=LTAI4FyeijYaddJVHKYZ2yg7&Signature=6KgeXSOn1sQAFN0nCcII8v8k2D4%3D" -O first_data.tar
tar -xvf first_data.tar
mv 20260105_ZhongGuoKeXueYuanBeiJingJiYinZuYanJiuSuoGuoJiaShengWuXinXiZhongXin-wubingqi-1_1 20260105_first
cd 20260105_first
md5sum -c MD5.txt
```

</details>




















