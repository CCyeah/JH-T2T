#################### 使用HiC将泛序列定位到染色体水平
### 
#①使用bowtie2-build为基因组建立索引

bowtie2-build  JH_t2t_6_chronly.fa JH_t2t_6_chronly

#②产生酶切位点信息.bed文件
python3 digest_genome.py -r mboi -o genome_mboi.bed JH_t2t_6_chronly.fa &
#③产生基因组大小信息文件.sizes文件
conda activate anno
seqkit fx2tab -l -n -i JH_t2t_6_chronly.fa  > genome.sizes
cat genome.sizes

#④修改HiCPro配置文件
cp ./*_1.fq.gz ./0_hic/fastq/ &
cp ./*_2.fq.gz ./0_hic/fastq/ &
gunzip *.fq.gz
nohup cat *_1.fq >./JHHIC_R1.fastq &
nohup cat *_2.fq >./JHHIC_R2.fastq &
gzip JHHIC_R1.fastq
gzip JHHIC_R2.fastq

#fastq 
nohup ./HiC-Pro  -i /disk202/caocy/5_JH_T2T/0_hic/fastq -o /disk202/caocy/5_JH_T2T/0_hic/out/  -c config-hicpro.txt &
 
nohup ./HiC-Pro \
-s proc_hic \
-i /disk202/caocy/5_JH_T2T/0_hic/fastq \
-o /disk202/caocy/5_JH_T2T/0_hic/out \
-c config-hicpro.txt &

nohup ./HiC-Pro \
-s proc_hic \
-i /disk202/caocy/5_JH_T2T/0_hic/out/bowtie_results/bwt2/ \
-o /disk202/caocy/5_JH_T2T/0_hic/out \
-c config-hicpro.txt &
 
 
 nohup ./HiC-Pro \
-s quality_checks \
-i /disk202/caocy/5_JH_T2T/0_hic/out/bowtie_results/bwt2/ \
-o /disk202/caocy/5_JH_T2T/0_hic/out \
-c config-hicpro.txt &
 
 nohup ./HiC-Pro \
-s merge_persample \
-i /disk202/caocy/5_JH_T2T/0_hic/out/hic_results/data/ \
-o /disk202/caocy/5_JH_T2T/0_hic/out \
-c config-hicpro.txt &
 
###合并
LANG=en; sort -T tmp -S 50% -k2,2V -k3,3n -k5,5V -k6,6n -m ./Sample1/Sample1.allValidPairs ./Sample2/Sample2.allValidPairs ./Sample3/Sample3.allValidPairs ./Sample4/Sample4.allValidPairs | awk -F"\t" 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=$2 || c2!=$5 || s1!=$3 || s2!=$6){print;c1=$2;c2=$5;s1=$3;s2=$6}' > /disk202/caocy/5_JH_T2T/0_hic/out/hic_results/data/sample_all/sample_all.allValidPairs

 ###########
  nohup ./HiC-Pro \
-s build_contact_maps \
-i /disk202/caocy/5_JH_T2T/0_hic/out/hic_results/data/ \
-o /disk202/caocy/5_JH_T2T/0_hic/out \
-c config-hicpro.txt &

  nohup ./HiC-Pro \
-s ice_norm   \
-i /disk202/caocy/5_JH_T2T/0_hic/out/hic_results/matrix/ \
-o /disk202/caocy/5_JH_T2T/0_hic/out \
-c config-hicpro.txt &

###########

python /disk202/caocy/software/HiCPlotter/HiCPlotter.py -f Sample_all_40000_iced.matrix -o hic_plot_40000 -r 40000 -tri 1 -bed sample_all_40000_abs.bed -n JH-t2t -wg 1 -chr chr20 -ext pdf &



