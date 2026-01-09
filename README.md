# Spatial-Transcriptomics
Learning some analysis code.

## 1.**SpaceRanger**
<details>
<summary> </summary>


</details>



## 2.**AnnData data structure**
<details>
<summary> </summary>

![AnnData](https://raw.githubusercontent.com/weiwei4396/Spatial-Transcriptomics/main/picture/anndata.jpg)
scverse 是一个专注于生命科学基础工具的组织和生态系统，最初聚焦于单细胞数据分析。它的优势在于出色的扩展性、灵活性以及与现有Python数据科学和机器学习工具的强大互操作性。

在scverse生态系统中, AnnData是用来将数据矩阵与这些注释关联起来的核心工具。为了提高效率, AnnData支持稀疏矩阵 (sparse matrices) 和部分读取(partial reading), 这样可以更快地处理大规模数据。AnnData在功能上与R生态系统中的数据结构 (比如Bioconductor的SummarizedExperiment或Seurat对象)相似, 但R包通常使用转置后的特征矩阵 (基因 x 细胞)。

在AnnData的核心中, 存储了一个稀疏或密集矩阵 (在scRNA-seq中就是计数矩阵), 称为X, 这个矩阵的维度是 obs_names x var_names (细胞 x 基因), 其中obs(观测值)对应细胞条形码, var(变量)对应基因标识符。 矩阵X被两个Pandas数据框(DataFrame)包围。其中obs保存细胞的注释信息, var保存基因的注释信息。 

AnnData还可以储存很多额外信息. 比如, 其他关于观测值和变量的多维数据(如UMAP) 储存在obsm和varm中, 图结构(比如细胞之间的关系或基因之间的关系)存储在obsp和varp中, 任何不适合其他槽位的非结构化数据都可以存储在uns中, 还可以通过layers存储矩阵X的额外值。例如, 可以在名为counts的层中存储未经标准化的原始计数数据, 而在默认的层中存储标准化后的数据。

[AnnData](https://mp.weixin.qq.com/s/0OFRSB3BZcNltHkp_1VG1Q)

- obs: 'in_tissue', 'array_row', 'array_col', 'gender', 'age', 'tissue', 'Strain', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'n_counts', 'n_genes'
- var: 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
- uns: 'spatial', 'log1p', 'hvg', 'pca', 'neighbors', 'umap', 'tsne', 'tissue_colors'
- obsm: 'spatial', 'X_pca', 'X_umap', 'X_tsne'
- varm: 'PCs'
- obsp: 'distances', 'connectivities'


</details>


## 3.**SpaceRanger outs folder**
<details>
<summary> </summary>

![outs](https://github.com/weiwei4396/Spatial-Transcriptomics/blob/main/picture/SpaceRanger_outs.jpg)

使用 [SpaceRanger](https://www.cnblogs.com/huanping/p/16839765.html)

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
BarcodeX+BarcodeY, 两种Barcode
<details>
<summary> </summary>

- 1.数据预处理, 提取Read1的barcode和UMI, 得到只由Barcode和UMI组成序列; 根据实验设计是否有barcodeZ, 提取不同的索引; 三宫格或九宫格: BarcodeX[1-8], BarcodeY+UMI[27-46]; 整个切片: BarcodeX[1-8], BarcodeY[27-34], BarcodeZ+UMI[53:72];

```shell
# 设置输入文件; sampath为当前测序的R1.fastq.gz和R2.fastq.gz所在的文件夹;
sampath=/data/database/MAGIC-seq-NG/Olfb
sample=OlfBulb
fastq1=${sampath}/${sample}_R1.fastq.gz
fastq2=${sampath}/${sample}_R2.fastq.gz

# whitelist为所有的barcode的组合, 如X和Y共4900种barcode, 每个barcode占一行;
whitelist=/data/database/MAGIC-seq-NG/Olfb/Mouse_Adult_Organ_T9_70_50um/whitelist
ID=/data/database/MAGIC-seq-NG/Olfb/Mouse_Adult_Organ_T9_70_50um/T9-ids-barcode.txt

# 建立的输出文件夹;
st_path=${sampath}/${sample}

# 参考基因组, 小鼠/人;
# 创建STAR的参考基因组index;
# /data/workdir/panw/software/STAR-2.7.11b/bin/Linux_x86_64/STAR --runThreadN 16 --runMode genomeGenerate --genomeDir star2710b --genomeFastaFiles /data/workdir/panw/reference/human/refdata-gex-GRCh38-2024-A/fasta/genome.fa --sjdbGTFfile /data/workdir/panw/reference/human/refdata-gex-GRCh38-2024-A/genes/genes.gtf --sjdbOverhang 100 --genomeSAindexNbases 14 --genomeChrBinNbits 18 --genomeSAsparseD 3
MAP=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/star2710b
ANN=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/genes/genes.gtf

log_file=${sampath}/${sample}_st_log.txt
t_num=32

# 注意: seqkit需要版本 2.0.0, 至少已知 2.10.0 中 concat功能不能支持现在的分析;
# 2.seqkit分割文件处理; 默认分成10份; 然后分别提取 BarcodeX和BarcodeY-UMI;
# 创建输出文件目录;
mkdir ${st_path}
mkdir ${st_path}/split
mkdir ${st_path}/out

FileNum=10
seqkit split2 -1 ${fastq1} -2 ${fastq2} -p $FileNum -O ${st_path}/split -f
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




# 使用自己写的脚本合并;
python concat_me.py -a test1-8.fastq.gz -b test27-46.fastq.gz -o BarcodeUMI_R1.fastq.gz

# 使用STARsolo将数据mapping;
/data/workdir/panw/software/STAR-2.7.11b/bin/Linux_x86_64/STAR --genomeDir ${MAP} \
  --outFileNamePrefix ${st_path}/STARsolo/${sample} \
  --readFilesCommand cat \
  --readFilesIn /data/database/MAGIC-seq-NG/Olfb/OlfBulb/Olf_R2.fastq /data/database/MAGIC-seq-NG/Olfb/OlfBulb/Olf_R1.fastq \
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

# 将STARsolo生成的基因表达转化为Anndata;




```


</details>


BarcodeX+BarcodeY+BarcodeZ, 三种Barcode
<details>
<summary> </summary>

</details>



## **Question & Answer**
<details>
<summary> </summary>

### 1. Segment Anything Model (SAM)



### 2. image registration (图像配准)

- 图像配准就是把两张不同来源, 不同大小, 不同角度的图像对齐, 使它们的同一个位置在像素空间上对应起来;
- H&E图是经过扫描仪拍摄得到的, 随着组织样本放置, 扫描角度、形变, 不同批次 H&E 图像中的网格, 不一定是正交的、等距的;
- 在H&E图和对应的spatial spot图上的同一个位置选取参考点, 两个对应点通常用于仿射变换, 线性缩放, 旋转, 平移;
- 像素坐标: 每个像素有自己的位置, 左上角永远是 (0, 0), 向右为x+, 向下为y+;






</details>




















