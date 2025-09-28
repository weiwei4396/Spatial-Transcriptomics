
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












