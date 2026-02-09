# Spatial-Transcriptomics
Learning some analysis code.

## 0.**Introduction**
### Sequencing-based spatial transcriptomics
<details>
<summary> </summary>

基于测序的空间组学平台通过将组织切片置于载玻片上，在一系列预设的空间测量位置捕获分子信息（这些信息可代表基因表达、DNA结合位点、抗体偶联标签等）。每个测量位置都带有唯一的空间条形码，用于标记其空间坐标；测序读段随后会根据基因或基因组区域等特征进行汇总（例如以计数形式表示）。

测序技术的主要优势在于，其所检测的分子特征通常是非靶向的（untargeted），即无需预先设计和优化探针面板。然而在实际应用中，许多空间组学实验仍需使用特定的探针组合——例如空间版CITE-seq；而像 Visium（v2 WT）和 Visium HD（WT）等平台虽号称“全转录组”，但其依赖于针对已知转录组设计的捕获探针面板，因此通常不适用于非模式生物。

不同平台的空间分辨率存在差异，主要取决于各测量位置的尺寸大小与间距。在具体生物样本中，由于组织细胞密度不同，每个空间测量位置可能包含零个、一个或多个细胞。这类平台中的测量单元通常被称为“spots”、“beads”或“bins”。

</details>

### Imaging-based spatial transcriptomics
<details>
<summary> </summary>

基于成像的平台（也称为基于单分子检测的平台）通过原位杂交（ISH）或原位测序（ISS）技术，对预设的靶向基因面板（通常包含数百至数千个基因）进行逐轮检测，从而在组织中原位识别单个RNA分子的空间位置。由于每个转录本都被单独识别，原始数据天然具有亚细胞级别的空间分辨率。

在预处理阶段，需通过图像分割（image segmentation）来识别单个细胞或细胞核的边界，并将检测到的 RNA分子分配给相应的细胞或细胞核。然而，细胞分割极具挑战性，尤其在组织切片中：由于细胞在三维空间中存在重叠，在二维切片的某一 (x, y) 位置可能包含多个细胞的信号，导致边界难以准确界定。

完成分割后，可将RNA分子按细胞汇总为细胞级别的基因表达计数，也可直接在单分子水平开展分析。若采用细胞级别数据，许多分析方法可复用自基于“spot”的空间转录组学或单细胞 RNA 测序领域已有的工具。

在实验设计阶段，一个关键环节是探针面板设计，即选择一组具有生物学意义的靶向基因。目前已有多种商业化靶向基因组合可供选择，适用于特定的生物学研究场景。

</details>


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


MAGIC-seq数据中对建库问题的排查
<details>
<summary> </summary>
0.通过每个位置碱基含量百分比, 初步判断read1和read2的结果;
```shell
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
```

1.按照确定位置找barcode, 查看正确的barcode比例; 通过改变-m参数可以查看m个barcode是否正确;
```shell
conda activate py310
```
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




## 5.**空间转录组分析模式**

#### 1.单变量分布模式：分析关注单个变量的空间分布, 而不考虑与其他变量的关系
- 意义在于识别并量化单个生物特征(如一个基因、一种细胞类型)在组织空间中的分布规律, 它解决了"它在哪里表达/存在? 分布是否均匀? 是否有空间趋势或聚集? "的问题, 这是后续所有复杂分析的基础;
- 空间趋势分析：使用广义加性模型或空间平滑方法, 模拟基因表达或细胞丰度随空间坐标(如从皮层表层到深层、从癌巢到间质)的连续变化梯度;
#### 2.二元空间关系
- 探究两个生物实体在空间上的相互作用关系。其核心是量化空间邻近性是否超过了随机预期, 从而推断潜在的功能性相互作用或空间排斥;
- 关键应用: 1.细胞-细胞相互作用: 量化免疫细胞与肿瘤细胞的共定位(如CD8+ T细胞浸润), 或细胞类型的相互排斥(如特定基质细胞与肿瘤区域分离); 2.配体-受体共定位: 评估配体表达细胞与受体表达细胞在空间上的邻近性, 为细胞通信提供空间证据; 3. 基因-基因共表达: 分析两个基因在单细胞或亚细胞水平的空间共表达模式, 提示潜在的共调控或功能协同;
- 实现方法: 1. 基于距离的共定位量化: 交叉K函数/Ripley's K：计算在给定半径范围内, 一种细胞类型(或基因阳性细胞)周围, 另一种细胞类型(或基因阳性细胞)的数目是否显著多于随机分布; A. 最近邻距离分析: 计算细胞类型A到细胞类型B的最近邻距离分布, 并与通过空间位置随机化生成的零分布进行比较; B. 细胞邻接分析: 基于Delaunay三角剖分或Voronoi图构建细胞邻接网络, 统计不同类型细胞互为邻居的频率;
- 实现方法: 2. 配体-受体共表达分析: 空间相关法：计算配体基因和受体基因表达向量的空间相关性(如Spearman相关); 细胞邻近通信评分: 对于每个细胞, 根据其邻近细胞中配体/受体的表达水平, 计算一个潜在的通信得分(如NicheNet、CellChat的空间扩展方法);
#### 3.高级结构
- 识别由多个变量(基因、细胞类型)协同定义的、在空间上连续或重复出现的功能单元或微环境。它超越了成对关系, 旨在发现组织中稳定且具有生物学意义的空间模块;
- 关键应用：定义空间功能域、细胞生态位、肿瘤微环境亚型、组织结构相关的多细胞群落;
- 实现方法：无监督空间聚类：基于特征的聚类：利用所有基因的表达谱(或重要基因子集)进行聚类(如K-means, 层次聚类), 并将聚类标签映射回空间, 检查其空间连续性; 常用工具如BayesSpace、SpaGCN专门优化了空间平滑性; 基于图的聚类：首先构建细胞的空间邻域图(基于坐标), 然后在图上进行聚类(如Leiden算法), 使聚类结果天然具有空间连续性;
- 实现方法：生态位/微环境识别：多元共现模式分析：利用非负矩阵分解或主题模型(如LDA), 将每个空间位置分解为若干"生态位主题"的混合, 每个主题由一组共现的细胞类型或基因定义; 基于细胞邻域组成的聚类：对于每个细胞, 统计其周围一定半径内所有细胞类型的组成, 形成一个"微环境向量", 然后对这些向量进行聚类, 从而定义重复出现的多细胞生态位;

参考

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




















