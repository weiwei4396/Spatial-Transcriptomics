
## 1.**SpaceRanger**
[官方网址](https://www.10xgenomics.com/support/software/space-ranger/latest)

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









