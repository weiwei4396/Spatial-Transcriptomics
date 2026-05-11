## decoder-seq脚本

必须准备的文件：
- barcode_coordinate.txt: 记录了barcode的顺序, 比如从X1到X75和从Y1到Y75, 以及对应的序列。必须包含三列"barcode"表示所有X和Y的组合序列(16bp), "x"表示x轴上排序的位置, 示例芯片是从1到75, "y"表示y轴上排序的位置, 示例芯片是从1到75。
- 全分辨率的H&E染色图。

分析流程主要包括三步：
- 1.根据barcode清单提取原始测序数据中有效的barcode。
- 2.将有效barcode的reads比对到参考基因组。
- 3.
首先





测试过这个分析流程，能根据编辑距离很好的校正barcode，缺点是相对较慢，与我自己写的脚本的区别在于它能校正1bp的插入和删除的barcode错误，我放弃了这部分，只保留了1bp的错配的barcode。而MAGIC-seq更粗暴，只是按位置提取出来，只保留了正确的barcode。


