# Spatial-Transcriptomics
Learning some analysis code.

1.AnnData data structure
<details>
<summary>Click to expand</summary>

![AnnData](https://raw.githubusercontent.com/weiwei4396/Spatial-Transcriptomics/main/picture/anndata.jpg)
AnnData是用来将数据矩阵与这些注释关联起来的核心工具。为了提高效率, AnnData支持稀疏矩阵 (sparse matrices) 和部分读取(partial reading), 这样可以更快地处理大规模数据。AnnData在功能上与R生态系统中的数据结构 (比如Bioconductor的SummarizedExperiment或Seurat对象)相似, 但R包通常使用转置后的特征矩阵 (基因 x 细胞)。

在AnnData的核心中, 存储了一个稀疏或密集矩阵 (在scRNA-seq中就是计数矩阵), 称为X, 这个矩阵的维度是 obs_names x var_names (细胞 x 基因), 其中obs(观测值)对应细胞条形码, var(变量)对应基因标识符。 矩阵X被两个Pandas数据框(DataFrame)包围。其中obs保存细胞的注释信息, var保存基因的注释信息。 

AnnData还可以储存很多额外信息. 比如, 其他关于观测值和变量的多维数据(如UMAP) 储存在obsm和varm中, 图结构(比如细胞之间的关系或基因之间的关系)存储在obsp和varp中, 任何不适合其他槽位的非结构化数据都可以存储在uns中, 还可以通过layers存储矩阵X的额外值。例如, 可以在名为counts的层中存储未经标准化的原始计数数据, 而在默认的层中存储标准化后的数据。
</details>


