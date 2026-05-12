
# 处理二代的MAGIC-seq流程; 以小鼠嗅球组织为例, T9芯片;
# ----------------------------------------------------------------------
# 1.准备数据;
# 文件夹和样品;
sampath=/data/database/MAGIC-seq-NG/Olfb/resultYihuan
sample=CRR1158889

fastq1=${sampath}/${sample}_R1_trim.fastq.gz
fastq2=${sampath}/${sample}_R2_trim.fastq.gz

# 给定的白名单; 
whitelist=/data/database/MAGIC-seq-NG/Olfb/Mouse_Adult_Organ_T9_70_50um/whitelist
#ID=/data/database/MAGIC-seq-NG/Olfb/Mouse_Adult_Organ_T9_70_50um/T9-ids-barcode.txt

MAP=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/star2710b
ANN=/data/workdir/panw/reference/mouse/refdata-gex-GRCm39-2024-A/genes/genes.gtf

# 线程;
t_num=16

ulimit -n 100000
mkdir ${sampath}/STARsolo

/data/workdir/panw/software/STAR-2.7.11b/bin/Linux_x86_64/STAR --genomeDir ${MAP} \
  --outFileNamePrefix ${sampath}/STARsolo/${sample}_ \
  --readFilesCommand zcat \
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

# zcat 可以换成 pigz -p 16 -dc











