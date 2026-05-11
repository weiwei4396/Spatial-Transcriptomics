## decoder-seq脚本

分析流程主要包括三步：
- 1.根据barcode清单提取有效的barcode。
- 2.将有效的reads比对到参考基因组。
- 3.


测试过这个分析流程，能根据编辑距离很好的校正barcode，缺点是相对较慢，与我自己写的脚本的区别在于它能校正1bp的插入和删除的barcode错误，我放弃了这部分，只保留了1bp的错配的barcode。而MAGIC-seq更粗暴，只是按位置提取出来，只保留了正确的barcode。


