import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import scanpy as sc
import warnings

# 设置绘图参数;
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["font.sans-serif"] = "Arial"
warnings.filterwarnings("ignore")
sc.settings.set_figure_params(dpi=100, dpi_save=300, facecolor="white")

# 1.读取文件;
# 载入Spaceranger处理visium后的文件夹, 里面包含spatial文件夹和h5文件;
h5_path = "/data/workdir/panw/Data/NAR_Mouse_spatial/MOB/GSM4656181_MOB-illumina"
# 读取h5文件, 参数count_file 默认为 "filtered_feature_bc_matrix.h5文件";
# 参数 source_image_path 默认为全分辨率的图片;
# 这个函数依次读取文件, 然后存储到AnnData中;
adata = sc.read_visium(h5_path)


# 2.首先对一些基因在组织上绘制表达;
# 设置图形的分辨率/大小, 样式和格式;
sc.set_figure_params(facecolor="white", figsize=(6,6), dpi_save=300, dpi=100)
# 在图像上叠加数据展示, img_key表示不同的图像, color绘制细胞或基因, cmap是色系;
# size是点的大小, bw表示图像是否是灰色图像, alpha_img表示绘制图像的alpha;
sc.pl.spatial(adata, img_key="hires", color=['in_tissue', 'Kctd12', 'Fabp7', 'Pcp4', 'Doc2g'], 
                cmap='Spectral_r', ncols=3, basis='spatial', size=1.5, 
                bw=False, alpha_img=1.)
# 绘制组织图像;
sc.pl.spatial(adata, img_key="hires", color=['in_tissue'], cmap='Spectral_r', spot_size=0)
# 关注切片中感兴趣的区域绘图; 聚焦图中的一小部分;
# crop_coord需要手动设置; adata.obsm['spatial']存储的是二位空间坐标信息;
# adata.obsm['spatial'].min(axis=0) 表示空间坐标最小的坐标, 在图的最左侧和最上侧;
# adata.obsm['spatial'].max(axis=0) 表示空间坐标最大的坐标, 在图的最右侧和最下测;
# crop_coord里面的四个元素是[左,右,上,下];
# vmax表示每个小图中色阶上限的值; 热图的点的颜色;
sc.pl.spatial(adata, img_key="hires",color=['Kctd12', 'Fabp7', 'Pcp4', 'Doc2g'], cmap='Spectral_r', ncols=3,
              basis='spatial', size=1.5, bw=False, alpha_img=1, crop_coord=[15000, 16000, 15000, 16000], 
              vmax=[15, 30, 6, 10])


# 3.探索数据结构;
# uns是字典数据类型, h5文件中的library id, 切片名称; 我这个例子是dict_keys(['164315']);
adata.uns['spatial'].keys()
# 每个切片下面存储的信息, 包括图片、缩放因子和metadata; dict_keys(['images', 'scalefactors', 'metadata'])
adata.uns['spatial']['164315'].keys()
# 图片就是包含高分辨率和低分辨率的图片; dict_keys(['hires', 'lowres'])
adata.uns['spatial']['164315']['images'].keys()
# 缩放因子存储的就是json文件中的内容, 包括spot半径, 基准点半径, 高低分辨率图片的scalefactors;
adata.uns['spatial']['164315']['scalefactors']
# metadata存储的是'chemistry_description': 'custom', 'software_version': '4509.7.5';
adata.uns['spatial']['164315']['metadata']


# 4.下面一直到聚类分析, 经典的单细胞数据分析流程;
# 首先展示对所有数据的总体展示;
# var中保存的是关于基因的注释信息; 计算线粒体基因的reads在细胞中的百分比;
adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
# 展示数据中的高表达基因;
sc.pl.highest_expr_genes(adata, n_top=20)

# 质控方面的指标; adata.var_names_make_unique()这一行是因为有些基因ID是有重复的, 这个函数将基因ID设置为不同;
adata.var_names_make_unique()
sc.set_figure_params(facecolor="white", figsize=(3,3), dpi_save=300, dpi=100)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, size=3, multi_panel=True)
# 查看质控指标之间的关系;
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8,4))
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show=False, ax=axes[0])
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=False, ax=axes[1])
plt.subplots_adjust(wspace=.3) # 用来调整子图之间的间距; wspace表示宽度方向的间距;


# 数据过滤;
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs['pct_counts_mt'] < 20]
print(f'#cells after MT filter:{adata.n_obs}')
sc.pp.filter_genes(adata, min_cells=10)

# 归一化;
# 对数据进行全局归一化，使得每个细胞的总表达值相等或接近;
sc.pp.normalize_total(adata, inplace=True)
# 对数据进行对数转换，以减小不同细胞之间的表达值范围差异;
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
# 展示高变基因;
sc.set_figure_params(facecolor="white", dpi_save=300, dpi=100)
sc.pl.highly_variable_genes(adata)
# 后续使用2000个高变基因;
adata.raw = adata
adata = adata[:,adata.var.highly_variable]

# 降维聚类;
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters")
# 绘制聚类图, 一个是spot的umap图, 一个是在图片上的聚类的图, 不考虑空间信息的聚类;
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12,6))
sc.pl.umap(adata, color=["clusters"], wspace=0.4, ax=axes[0], show=False)
sc.pl.spatial(adata, img_key="hires", color=["clusters"], ax=axes[1], show=False)
plt.subplots_adjust(wspace=0.3)
# 绘制经典marker在umap聚类图上的展示;
sc.pl.umap(adata, color=['Sulf1'], cmap='Spectral_r', ncols=3)


