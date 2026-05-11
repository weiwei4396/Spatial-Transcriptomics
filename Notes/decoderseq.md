### 4.4 decoder-seq的脚本
[Decoder-seq分析流程](https://github.com/songjiajia2018/Decoder-seq/tree/main)

测试过这个分析流程，能根据编辑距离很好的校正barcode，缺点是相对较慢，与我自己写的脚本的区别在于它能校正1bp的插入和删除的barcode错误，我放弃了这部分，只保留了1bp的错配的barcode。而MAGIC-seq更粗暴，只是按位置提取出来，只保留了正确的barcode。


