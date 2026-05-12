## Content
- [AnnData](#AnnData)
- [STARsolo](#STARsolo)
- [fastp](#fastp)

## 1.**AnnData**

[AnnData结构](https://raw.githubusercontent.com/weiwei4396/Spatial-Transcriptomics/main/picture/anndata.jpg)

scverse 是一个专注于生命科学基础工具的组织和生态系统，最初聚焦于单细胞数据分析。它的优势在于出色的扩展性、灵活性以及与现有Python数据科学和机器学习工具的强大互操作性。

在scverse生态系统中，AnnData是用来将数据矩阵与这些注释关联起来的核心工具。为了提高效率，AnnData支持稀疏矩阵 (sparse matrices) 和部分读取(partial reading)，这样可以更快地处理大规模数据。AnnData在功能上与R生态系统中的数据结构 (比如Bioconductor的SummarizedExperiment或Seurat对象)相似，但R包通常使用转置后的特征矩阵 (基因 x 细胞)。

在AnnData的核心中，存储了一个稀疏或密集矩阵 (在scRNA-seq中就是计数矩阵)，称为X，这个矩阵的维度是 obs_names x var_names (细胞 x 基因)，其中obs(观测值)对应细胞条形码，var(变量)对应基因标识符。 矩阵X被两个Pandas数据框(DataFrame)包围。其中obs保存细胞的注释信息，var保存基因的注释信息。 

AnnData还可以储存很多额外信息。比如，其他关于观测值和变量的多维数据(如UMAP) 储存在obsm和varm中，图结构(比如细胞之间的关系或基因之间的关系)存储在obsp和varp中，任何不适合其他槽位的非结构化数据都可以存储在uns中，还可以通过layers存储矩阵X的额外值。例如，可以在名为counts的层中存储未经标准化的原始计数数据，而在默认的层中存储标准化后的数据。

[AnnData参考](https://mp.weixin.qq.com/s/0OFRSB3BZcNltHkp_1VG1Q)

<details>
<summary> AnnData数据存储信息结构 </summary>

- obs: 'in_tissue', 'array_row', 'array_col', 'gender', 'age', 'tissue', 'Strain', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'n_counts', 'n_genes'
- var: 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
- uns: 'spatial', 'log1p', 'hvg', 'pca', 'neighbors', 'umap', 'tsne', 'tissue_colors'
- obsm: 'spatial', 'X_pca', 'X_umap', 'X_tsne'
- varm: 'PCs'
- obsp: 'distances', 'connectivities'

</details>


## 2.**STARsolo**
[STARsolo官方教程](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)

- STARsolo的基因组索引与正常STAR运行相同。
- STARsolo中的read1完全是barcode+umi，且在输入参数中read2.fastq.gz放在前面，read1.fastq.gz放在后面。
- 为了使结果与cellranger相同，可以使用cellranger的注释，生成参考基因组索引。
- 细胞过滤参数--soloCellFilter，EmptyDrops_CR是[EmptyDrop算法](https://link.springer.com/article/10.1186/s13059-019-1662-y)
- 不同转录组特征的量化可以使用参数--soloFeatures，Gene是落在外显子区域的所有reads的计数，GeneFull是pre-mRNA计数，统计所有与基因位点重叠的reads，包括外显子和内含子reads，一般适合单核RNA-seq的分析。Velocyto计算每个细胞中每个基因的剪接、未剪接和模糊基因计数。SJ表示已注释和新发现的剪接连接点的计数。
- 多基因读取参数--soloMultiMappers，有多种算法可以使用，我常用的是EM算法。
- BAM标签是参数--outSAMattributes，CR/UR表示未校正的CellBarcode/UMI序列，CY/UY表示CellBarcode/UMI的质量评分，GX/GN表示基因ID和symbol ID，CB/UB表示校正后的CellBarcode/UMI。
- 对UMI去重使用--soloUMIdedup参数，有多种算法可以使用，1MM_All算法是把彼此之间只差1个碱基的UMI合并。1MM_Directional_UMItools，不只是看UMI是否差1个碱基，还会考虑每个 UMI的read支持数。1MM_Directional方法与前者类似，但更保守。Exact只合并完全相同的UMI，NoDedup不做UMI去重，1MM_CR模拟Cell Ranger 2–4版本的1 mismatch UMI collapsing算法。我现在用的是1MM_All。
```shell
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

























