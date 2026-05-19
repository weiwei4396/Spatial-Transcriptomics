

import script3_magicseq_downstream as magic
import os
import warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use('Agg')

def main_barcode2():
    # 嗅球样品;
    sample='sample1-7'
    chip_type='T9' # 芯片类型;
    reg='reg7' # 区域
    channels_num=70 # Number of channels
    barcode_numA=70 # Number of barcodes
    barcode_numB=70 # Number of barcodes
    res_um=50 # resolution

    Barcode_file_path='/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/data/Mouse_Adult_Organ_T9_70_50um/' # 所有barcode信息的文件夹;
    file_path='/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/Olfb/' # STARsolo比对后的文件夹;
    image_file_path='/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/data/T9-sample1-7/' # 所有图片信息的文件夹;
    geneinfo_path='/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/data/geneInfo.txt' # 基因信息文件;
    vit_path = '/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/data/sam_vit_h_4b8939.pth'

    output_path = os.path.join(file_path, 'pictures')
    os.makedirs(output_path, exist_ok=True)

    # Point coordinates for image registration
    HE_point1=(420,328) # HE1
    Spot_point1=(408,234) # marker1
    HE_point2=(4688,4653) # HE2
    Spot_point2=(4694,4582) # marker2

    line_point_r1c1=(713,704)  # spot1
    line_point_r1c70=(741,4221) # spot2
    line_point_r70c1=(4231,660) # spot3

    adata = magic.get_adata_STARsolo(sample=sample, chip_type=chip_type, reg=reg, 
                                    channels_num=channels_num, barcode_numA=barcode_numA, barcode_numB=barcode_numB,
                                    res_um=res_um, EM=False, Velocyto=False, soloFeatures='Gene',
                                    raw=False, species='mouse', 
                                    file_path=file_path, image_file_path=image_file_path, 
                                    Barcode_file_path=Barcode_file_path, geneinfo_path=geneinfo_path, 
                                    sam_checkpoint=vit_path,
                                    HE_point1=HE_point1, Spot_point1=Spot_point1,
                                    HE_point2=HE_point2, Spot_point2=Spot_point2,
                                    line_point_r1c1=line_point_r1c1, line_point_r1c70=line_point_r1c70, line_point_r70c1=line_point_r70c1)

    adata_filtered = magic.adata_filter_norm_process(adata, outputPath=output_path, cluster_key='louvain')
    adata_filtered.write_h5ad(f'{file_path}/{sample}_{chip_type}_{reg}.h5ad')


def main_barcode3():
    # E17
    sample='E17-1'
    chip_type='M9' # 芯片类型;
    reg='reg1' # 区域
    channels_num=450 # Number of channels
    barcode_numA=150 # Number of barcodes
    barcode_numB=150 # Number of barcodes
    res_um=20 # resolution

    Barcode_file_path='/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/E17-1/Barcode-M9-150-E17.5/' # 所有barcode信息的文件夹;
    file_path='/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/E17-1/' # STARsolo比对后的文件夹;
    image_file_path='/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/E17-1/M9-E17-1/' # 所有图片信息的文件夹;
    geneinfo_path='/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/data/geneInfo.txt' # 基因信息文件;
    vit_path = '/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/data/sam_vit_h_4b8939.pth'

    output_path = os.path.join(file_path, 'pictures')
    os.makedirs(output_path, exist_ok=True)

    # Point coordinates for image registration
    HE_point1=(97,178) # HE1
    Spot_point1=(162,242) # marker1
    HE_point2=(4708,4939) # HE2
    Spot_point2=(4817,4742) # marker2

    line_point_r1c1=(377, 311)  # spot1
    line_point_r1c70=(397, 4645) # spot2
    line_point_r70c1=(4574, 322) # spot3

    # 注意add_M9_15和add_M9_20参数, 可以填充缝隙;
    adata = magic.get_adata_STARsolo(sample=sample, chip_type=chip_type, reg=reg, 
                                    channels_num=channels_num, barcode_numA=barcode_numA, barcode_numB=barcode_numB,
                                    res_um=res_um, EM=False, Velocyto=False, soloFeatures='Gene',
                                    raw=False, species='mouse', 
                                    file_path=file_path, image_file_path=image_file_path, 
                                    Barcode_file_path=Barcode_file_path, geneinfo_path=geneinfo_path, 
                                    sam_checkpoint=vit_path,
                                    add_M9_15=False, add_M9_20=True,
                                    HE_point1=HE_point1, Spot_point1=Spot_point1,
                                    HE_point2=HE_point2, Spot_point2=Spot_point2,
                                    line_point_r1c1=line_point_r1c1, line_point_r1c70=line_point_r1c70, line_point_r70c1=line_point_r70c1)

    adata_filtered = magic.adata_filter_norm_process(adata, outputPath=output_path, cluster_key='louvain', plot_marker=False)
    adata_filtered.write_h5ad(f'{file_path}/{sample}_{chip_type}_{reg}.h5ad')


if __name__ == "__main__":
    main_barcode2()
    # main_barcode3()






