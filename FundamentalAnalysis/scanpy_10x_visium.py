"""
AnnData是用来将数据矩阵与这些注释关联起来的核心工具。
为了提高效率, AnnData支持稀疏矩阵 (sparse matrices) 和部分读取(partial reading), 这样可以更快地处理大规模数据.
AnnData在功能上与R生态系统中的数据结构 (比如Bioconductor的SummarizedExperiment或Seurat对象)相似,
但R包通常使用转置后的特征矩阵 (基因 x 细胞).

在AnnData的核心中, 存储了一个稀疏或密集矩阵 (在scRNA-seq中就是计数矩阵), 称为X, 这个矩阵的维度是 obs_names x var_names (细胞 x 基因),
其中obs(观测值)对应细胞条形码, var(变量)对应基因标识符. 矩阵X被两个Pandas数据框(DataFrame)包围.
其中obs保存细胞的注释信息, var保存基因的注释信息. 
AnnData还可以储存很多额外信息. 比如, 其他关于观测值和变量的多维数据(如UMAP) 储存在obsm和varm中, 
图结构(比如细胞之间的关系或基因之间的关系)存储在obsp和varp中,
任何不适合其他槽位的非结构化数据都可以存储在uns中, 还可以通过layers存储矩阵X的额外值.
例如, 可以在名为counts的层中存储未经标准化的原始计数数据, 而在默认的层中存储标准化后的数据。
"""

import scanpy as sc

# 载入Spaceranger处理visium后的文件夹, 里面包含spatial文件夹和h5文件;
h5_path = "/data/workdir/panw/Data/NAR_Mouse_spatial/MOB/GSM4656181_MOB-illumina"
# 读取h5文件, 参数count_file 默认为 "filtered_feature_bc_matrix.h5文件";
adata = sc.read_visium(h5_path)

# 首先对一些基因在组织上绘制表达;
# 设置图形的分辨率/大小, 样式和格式;
sc.set_figure_params(facecolor="white", figsize=(6,6), dpi_save=300, dpi=100)
# 在图像上叠加数据展示, img_key表示不同的图像, color绘制细胞或基因, cmap是色系;
# size是点的大小, bw表示图像是否是灰色图像, alpha_img表示绘制图像的alpha;
sc.pl.spatial(adata, img_key="hires", color=['in_tissue', 'Kctd12', 'Fabp7', 'Pcp4', 'Doc2g'], 
                cmap='Spectral_r', ncols=3, basis='spatial', size=1.5, 
                bw=False, alpha_img=1.)
# 绘制组织图像;
sc.pl.spatial(adata, img_key="hires", color=['in_tissue'], cmap='Spectral_r', spot_size=0)












