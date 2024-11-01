############################################################################################################
##基因组组装评估
# 统计组装基本结果
python ./contigStas_cc.py JH_t2t.fa  JH_t2t_contigStas.txt
############################################################################################################

bioawk -c fastx '{print ">"$name;print length($seq)}'  ./JH_t2t.fa > JH.summary
bioawk -c fastx '{print ">"$name;print length($seq)}'  ./Sscrofa11.fa > ref.summary

############################################################################################################
#busco
for j in `ls ./*.fasta`
do
prefix=`echo ${j} | awk -F "." '{print $j}'`
prefixid=`echo ${j} | awk -F "." '{print $j}'`
echo "
nohup busco -i ${prefix} -c 20 -o ${prefix} -m geno -l ./mammalia_odb10/ &
" >> busco.sh
done
nohup bash busco.sh &> busco.log &

############################################################################################################
#QUAST
quast.py ./JH_t2t.fa
quast.py ./Sus_scrofa.fasta
quast.py ./NX.fa
quast.py ./MS.fa

########### mergury QV
#2,694.52Mbp v 2679795048
$MERQURY/best_k.sh 2694520000
$MERQURY/best_k.sh 2694.52Mbp
#genome: 2694520000000
#tolerable collision rate: 0.001
#25.6288
#1. Build meryl dbs
    #meryl k=25 count output read1.meryl /disk191_3/LabData/Pig/1_DNAseq/Third_1/0_illumina/illumina/clean/FDES220007339-1a_L1_1_clean.rd.fq.gz
	nohup meryl k=20 count output read1.meryl./FDES220007339-1a_L1_1_clean.rd.fq.gz &
	nohup meryl k=20 count output read2.meryl./FDES220007339-1a_L1_2_clean.rd.fq.gz &
	nohup meryl k=20 count output read3.meryl./FDES220007339-1a_L2_1_clean.rd.fq.gz &
	nohup meryl k=20 count output read4.meryl./FDES220007339-1a_L2_2_clean.rd.fq.gz &
#Merge
meryl union-sum output child.meryl read*.meryl
#2. Build maternal meryl dbs 
mkdir maternal_sp
nohup meryl k=20 count output read1.meryl ./jin1/jin1_1.clean.fq.gz &
nohup meryl k=20 count output read2.meryl ./jin1/jin1_2.clean.fq.gz &
#Merge
meryl union-sum output maternal.meryl read*.meryl
#3. Build paternal meryl dbs 
mkdir paternal_sp
nohup meryl k=20 count output read1.meryl ./jin2/jin2_1.clean.fq.gz &
nohup meryl k=20 count output read2.meryl ./jin2/jin2_2.clean.fq.gz &
#Merge
meryl union-sum output paternal.meryl read*.meryl
$MERQURY/merqury.sh ./child.meryl ./JH_t2t.fa out_prefix_chr
#trio
$MERQURY/merqury.sh ./child.meryl ./paternal.meryl ./maternal.meryl ./JH_t2t.fa out_prefix_chr &
######

#Synteny analysis
############################################## mummer 共线性分析
#用只包含20条染色体的数据进行比较
 $reference 
 $query
nohup nucmer -t 36 -l 100 -c 1000 --prefix=JH $reference $query &
# 结果过滤 & 转换格式
for breed in JH
do
delta-filter -1 ${breed}.delta > ${breed}.filter
show-coords -rcl ${breed}.filter > ${breed}.coords
done

show-coords -TH ${breed}.filter | awk '$5>50000, OFS="\t"{print $8,$1,$2,$9,$3,$4}'> links_${breed}.tsv

nohup GenomeSyn -g1 JH_t2t_6_chronly.fa  -g2 Sscrofa11_1.fa \
-cen1 cen_t2t.bed -cen2 cen_sus.bed \
-cf1 ReferencevsQuery1.delta.filter.coords \
-GD1 T2T.gtf -GD2 SUS.gtf \
-PAV1 SV_t2T.bed -PAV2 SV_SUS11.1.bed &

#################################################################################
##snp  indel

echo "
minimap2 -cx asm5 -t8 --cs /disk222/caocy/1_T2T/7_ASS_2/0_basic/t2t_fa /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/Ssc.fa > asm.paf   &
sort -k6,6 -k8,8n asm.paf > asm.srt.paf           # sort by reference start coordinate
k8 paftools.js call asm.srt.paf -f /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa -s D > out.vcf
paftools.js call asm.srt.paf -f /disk222/caocy/1_T2T/5_gapfill_TGS_2/6_chr/V2/JH_t2t.fa -s D -L 20000 > out.vcf" >> minimap_GvG.sh
nohup bash minimap_GvG.sh &> minimap_GvG.log &


##########序列一致性评估###############################################
#使用比对工具minimap2（默认参数）比对回组装好的基因组
source /disk212/caocy/software/anaconda3/bin/activate /disk212/caocy/software/anaconda3/envs/python3.7
#hifi 
echo "
minimap2 -ax map-pb ${T2T} ${hifi_reads} -t 30 > aln_t2t.sam
samtools view -bS aln_t2t.sam > aln_t2t.bam
samtools sort -@ 20 aln_t2t.bam -O minimap_t2t_sorted.merged.bam --output-fmt BAM" >> minimap_t2t_hifi.sh

nohup bash minimap_t2t_hifi.sh &> minimap_t2t_hifi.log &
samtools sort -@ 30 -o minimap_t2t_sorted.merged.bam aln_t2t.bam
samtools index minimap_t2t_sorted.merged.bam

#ont 
echo "
minimap2 -ax map-ont ${T2T} ${ont_reads} -t 36 > aln_t2t.sam
samtools view -bS aln_t2t.sam > aln_t2t.bam
samtools sort -@ 30 aln_t2t.bam -O minimap_t2t_sorted.merged.bam --output-fmt BAM" >> minimap_t2t_ont.sh

nohup bash minimap_t2t_ont.sh &> minimap_t2t_ont.log &
samtools sort -@ 30 -o minimap_t2t_sorted.merged.bam aln_t2t.bam &
samtools index minimap_t2t_sorted.merged.bam &

samtools flagstat minimap.merged.bam > minimap.merged.bam.flagstat
samtools depth -aa minimap_t2t_sorted.merged.bam > depth.info


######################################################################################################

for i in fasqname
do
prefix=${i}
echo "
bwa mem -t 30 ref ./${prefix}.R1.fq.gz./${prefix}.R2.fq.gz > ${prefix}.sam
samtools view -b ${prefix}.sam > ${prefix}.bam
samtools sort -@ 30 -o ${prefix}.sorted.bam ${prefix}.bam
samtools index ${prefix}.sorted.bam
rm ${prefix}.sam ${prefix}.bam
" >> map.sh
done
nohup bash map_JH.sh &> map_JH.log &

for j in `ls *bam`
do
prefix=`echo ${j} | awk -F "." '{print $j}'`
samtools stats ${prefix} > ${prefix}.stats &
done

#Synteny analysis
############################################## mummer 共线性分析
### server 212
#conda create -n mummer4 -c conda-forge -c bioconda mummer4
#conda create -n mummer4
#conda install -c bioconda mummer4
#conda install gnuplot=4.6.0
source /disk212/caocy/software/anaconda3/bin/activate
conda activate mummer4
#用只包含20条染色体的数据进行比较
cat /disk212/miaoj/assemblies/rawData/Sscrofa11.fa | grep ">" | wc -l
reference=/disk212/miaoj/assemblies/rawData/Sscrofa11.fa
query=/disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa
#nohup nucmer -t 48 -l 100 -c 1000 --prefix=JH $reference $query &
#query=/disk212/miaoj/assemblies/rawData/JH.fa
nohup nucmer -t 36 -l 100 -c 1000 --prefix=JH $reference $query &
# 结果过滤 & 转换格式
for breed in JH
do
delta-filter -1 ${breed}.delta > ${breed}.filter
show-coords -rcl ${breed}.filter > ${breed}.coords
done

show-coords -TH ${breed}.filter | awk '$5>50000, OFS="\t"{print $8,$1,$2,$9,$3,$4}'> links_${breed}.tsv

mummerplot -p JH -f JH.filter 

mummerplot --png  -p test JH_hic.delta -R /disk212/miaoj/assemblies/rawData/Sscrofa11.fa  \
-Q /disk191_3/LabData/Pig/1_DNAseq/Third_1/1_ASS/X101SC22081679-Z01-F005/genome.review4.assembly.FINAL.fasta 
#A.genome.Chr.fa 只有染色体序列 (我最近学习的展示方式)

wget https://github.com/tpoorten/dotPlotly/blob/master/mummerCoordsDotPlotly.R
./mummerCoordsDotPlotly.R -i JH.coords -o out -s -t -m 500 -q 500000 -k 7 -l


 nucmer --mum -D 5 -t 16 $reference $query -p prefix
 delta-filter –i 90 prefix.delta -1 > prefix.best.delta
 mummerplot -p prefix -f prefix.best.delta -t postscript
 /usr/bin/ps2pdf prefix.ps prefix.pdf


source /disk212/caocy/software/anaconda3/bin/activate
conda activate mummer4
#git https://github.com/JM-SONG/GenomeSyn.git
cd /disk212/caocy/software/GenomeSyn/GenomeSyn-1.2.7/
source ./install.sh

##方便画图将染色名称由chr01改为1
cp /disk222/caocy/1_T2T/5_gapfill_TGS/17_JLcontig/result.scaff_seqs .
perl /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/split.pl result.scaff_seqs

perl /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/split.pl /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/Sscrofa11.fa
cat chr*.fa >Ssc.fa
cat JHt2t.fa | grep ">" | wc -l
cd /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/JHt2t/
bwa index -p index/JH_t2t /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/JHt2t.fa &

#sed -i "s#Chr02#2#g" chr2.fa
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 X Y
do
sed -i "s#chr$j#$j#g" chr$j.fa
done
cat chr*.fa >Ssc_n.fa

cd /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/JHt2t/
bwa index -p index/JH_t2t /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/JHt2t.fa &

nohup GenomeSyn -g1 /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa  -g2 /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/Ssc.fa \
-GD1 /disk222/caocy/1_T2T/6_Anno/0_lift/t2t.gene.gff3 -GD2 SUS.gff3 \
source /disk212/caocy/software/anaconda3/bin/activate /disk212/caocy/software/anaconda3/envs/T2T

GenomeSyn -g1 /disk212/miaoj/assemblies/rawData/Sscrofa11.fa  -g2 /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa
###################mannal
GenomeSyn -t 3 -n3 12 -g1 ../data/rice_MH63.fa -g2 ../data/rice_ZS97.fa -g3 ../data/rice_R498.fasta \
-cf1 ../data/rice_MH63vsZS97.delta.filter.coords -cf2 ../data/rice_MH63vsR498.delta.filter.coords \
#着丝粒
-cen1 ../data/rice_MH63_centromere.bed -cen2 ../data/rice_ZS97_centromere.bed -cen3 ../data/rice_R498_centromere.bed \  
#端粒 
-tel1 ../data/rice_MH63_telomere.bed -tel2 ../data/rice_ZS97_telomere.bed -tel3 ../data/rice_R498_telomere.bed \
#te
-TE2 ../data/rice_ZS97_repeat.bed -PAV1 ../data/rice_MH63_PAV.bed -PAV2 ../data/rice_ZS97_PAV.bed \
-NLR1 ../data/rice_MH63_NLR.bed -NLR2 ../data/rice_ZS97_NLR.bed \
-r MH63 -q1 ZS97 -q2 R498 \
-GD1 ../data/rice_MH63_nonTEgene.gff3 -GD2 ../data/rice_ZS97_nonTEgene.gff3 -GD3 ../data/rice_R498_IGDBv3_coreset.gff \
-GC2 ../data/rice_ZS97_GC_10000.bed -GC_win 100000 -TE_min 40


GenomeSyn -g1 /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa  -g2 /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/Ssc.fa


#wget https://ftp.ensembl.org/pub/release-109/gff3/sus_scrofa/Sus_scrofa.Sscrofa11.1.109.chr.gff3.gz
cd /disk222/caocy/1_T2T/7_ASS/3_mummer/GenomeSyn_3/


source /disk212/caocy/software/anaconda3/bin/activate /disk212/caocy/software/anaconda3/envs/T2T
#conda install -c bioconda bioawk

nohup GenomeSyn -g1 /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa  -g2 /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/Ssc.fa \
-GD1 /disk222/caocy/1_T2T/6_Anno/0_lift/t2t.gene.gff3 -GD2 SUS.gff3 \
-cf1 /disk222/caocy/1_T2T/7_ASS/3_mummer/GenomeSyn_2/JHvsSSC.delta.filter.coords&

cd /disk222/caocy/1_T2T/7_ASS/3_mummer/
nohup GenomeSyn -g1 /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa  \
-g2 /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/Ssc.fa \
-GD1 /disk222/caocy/1_T2T/6_Anno/0_lift/t2t.gene.gff3 \
-GD2 /disk222/caocy/1_T2T/6_Anno_2/Sus_scrofa.Sscrofa11.1.109.chr.gff3 \
-cf1 /disk222/caocy/1_T2T/7_ASS/3_mummer/GenomeSyn_2/JHvsSSC.delta.filter.coords&





########### 4_mergury###############基于二代测序、kmer计算QV值
/disk212/caocy/5_T2T/3_verkko/
#2,694.52Mbp v 2679795048
$MERQURY/best_k.sh 2694520000
$MERQURY/best_k.sh 2694.52Mbp
#genome: 2694520000000
#tolerable collision rate: 0.001
#25.6288
#1. Build meryl dbs
    #meryl k=25 count output read1.meryl /disk191_3/LabData/Pig/1_DNAseq/Third_1/0_illumina/illumina/clean/FDES220007339-1a_L1_1_clean.rd.fq.gz
	nohup meryl k=20 count output read1.meryl /disk191_3/LabData/Pig/1_DNAseq/Third_1/0_illumina/illumina/clean/FDES220007339-1a_L1_1_clean.rd.fq.gz &
	nohup meryl k=20 count output read2.meryl /disk191_3/LabData/Pig/1_DNAseq/Third_1/0_illumina/illumina/clean/FDES220007339-1a_L1_2_clean.rd.fq.gz &
	nohup meryl k=20 count output read3.meryl /disk191_3/LabData/Pig/1_DNAseq/Third_1/0_illumina/illumina/clean/FDES220007339-1a_L2_1_clean.rd.fq.gz &
	nohup meryl k=20 count output read4.meryl /disk191_3/LabData/Pig/1_DNAseq/Third_1/0_illumina/illumina/clean/FDES220007339-1a_L2_2_clean.rd.fq.gz &
#Merge
meryl union-sum output child.meryl read*.meryl
#2. Build maternal meryl dbs 母 222 213
mkdir maternal_sp
nohup meryl k=20 count output read1.meryl /disk191_3/LabData/Pig/1_DNAseq/Third_1/0_illumina/X101SC22073411-Z01-J019/00.CleanData/jin1/jin1_1.clean.fq.gz &
nohup meryl k=20 count output read2.meryl /disk191_3/LabData/Pig/1_DNAseq/Third_1/0_illumina/X101SC22073411-Z01-J019/00.CleanData/jin1/jin1_2.clean.fq.gz &
#Merge
meryl union-sum output maternal.meryl read*.meryl
#3. Build paternal meryl dbs 父202 201
mkdir paternal_sp
nohup meryl k=20 count output read1.meryl /disk191_3/LabData/Pig/1_DNAseq/Third_1/0_illumina/X101SC22073411-Z01-J019/00.CleanData/jin2/jin2_1.clean.fq.gz &
nohup meryl k=20 count output read2.meryl /disk191_3/LabData/Pig/1_DNAseq/Third_1/0_illumina/X101SC22073411-Z01-J019/00.CleanData/jin2/jin2_2.clean.fq.gz &
#Merge
meryl union-sum output paternal.meryl read*.meryl
$MERQURY/merqury.sh /disk212/caocy/5_T2T/3_verkko/child.meryl /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa out_prefix_chr
#trio
$MERQURY/merqury.sh /disk212/caocy/5_T2T/3_verkko/child.meryl /disk212/caocy/5_T2T/3_verkko/paternal.meryl /disk212/caocy/5_T2T/3_verkko/maternal.meryl /disk228/wz/caocy/1_T2T/5_polish/sgs/01_rundir/genome.nextpolish.fasta out_prefix_chr &

########### circos画图准备
### 染色体长度
#计算基因密度
source /disk212/caocy/software/anaconda3/bin/activate /disk212/caocy/software/anaconda3/envs/T2T
python /disk222/caocy/1_T2T/7_ASS/3_mummer/gc_gene_D.py





############################
#SNP
#reference=/disk212/miaoj/assemblies/rawData/Sscrofa11.fa
#query=/disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa
#ref.fa=/disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa
#asm.fa=/disk212/miaoj/assemblies/rawData/Sscrofa11.fa
source /disk212/caocy/software/anaconda3/bin/activate
conda activate tgsgapcloser 
nohup minimap2 -cx asm5 -t8 --cs /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa\
 /disk212/miaoj/assemblies/rawData/Sscrofa11.fa > asm.paf  &  # keeping this file is recommended; --cs required!
sort -k6,6 -k8,8n asm.paf > asm.srt.paf           # sort by reference start coordinate
k8 /disk212/caocy/software/anaconda3/envs/tgsgapcloser/bin/paftools.js call asm.srt.paf > asm.var.txt
k8 /disk212/caocy/software/anaconda3/envs/tgsgapcloser/bin/paftools.js call asm.srt.paf -f /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa -s D > out.vcf
paftools.js call asm.srt.paf -f /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa -s D -L 20000 > out.vcf

samtools faidx /disk191/chenzt/database/fasta/genome.fa

bedtools makewindows -g genome.fa.fai -w 1000000 >region.bed

bedtools coverage -a /disk222/caocy/1_T2T/7_ASS/1_mapping_hifi/50k.bed -b SNPposition.txt  -counts >SNP-density.txt
#minimap2 -ax asm5 /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa /disk212/miaoj/assemblies/rawData/Sscrofa11.fa > genome1_genome2.sam
#minimap2 -x asm5 genome1.fa genome2.fa > alignment.sam
#cd /disk222/caocy/software/
#git clone https://github.com/lh3/htsbox
#(cd htsbox && make)
minimap2 -axasm5 /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa\
 /disk212/miaoj/assemblies/rawData/Sscrofa11.fa | samtools sort - > sorted.bam
#htsbox/htsbox pileup -q5 -S10000 -vcf wt_minion.fasta sorted.bam > diff.vcf
#minimap2 -x splice /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa /disk212/miaoj/assemblies/rawData/Sscrofa11.fa > aln_t2t.sam
#minimap2 -x splice /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa /disk212/miaoj/assemblies/rawData/Sscrofa11.fa > aln_t2t.sam
#-axasm5
#samtools view -bS aln_t2t.sam > aln_t2t.bam
#samtools sort -@ 20 aln_t2t.bam -O minimap_t2t_sorted.merged.bam --output-fmt BAM
#/disk222/caocy/software/htsbox/htsbox pileup -q5 -S10000 -vcf /disk222/caocy/1_T2T/7_ASS/1_mapping_ngs/index/JH_t2t.fa aln_t2t.bam > diff.vcf
