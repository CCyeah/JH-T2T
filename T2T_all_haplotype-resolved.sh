##hifiasm_trio
source /disk212/caocy/software/anaconda3/bin/activate  T2T
# trio-binning模式，需要额外安装yak
conda install -c bioconda yak

cd 0_hifiasm_trio
ln -s ./jin*/jin*_*.clean.fq.gz .
#母
yak count -k 20  -b 37 -t 20 -o  IMCGF.yak <(zcat jin1_1.clean.fq.gz ) <(zcat jin1_2.clean.fq.gz) 
#父
yak count -k 20  -b 37 -t 20 -o  IMCGM.yak <(zcat jin2_1.clean.fq.gz ) <(zcat jin2_2.clean.fq.gz) 
nohup hifiasm -o JH -t 30 -1 IMCGF.yak -2 IMCGM.yak ./JH_hifi.fastq.gz > JH_hifi.log 2>&1 &



zcat ./onecellall.fastq.gz ./pass.fq.gz >ont.pass.fastq

nohup hifiasm -o JH -t 48 --ul ont.pass.fastq.gz  -1 ./0_hifiasm_trio/IMCGF.yak -2 ./0_hifiasm_trio/IMCGM.yak  --dual-scaf --telo-m CCCTAA ./0_hifi/5cell/JH_hifi.fastq.gz > JH_hifi.log 2>&1 &

#nohup bam2fastq -o JH_ONT ./1_ASS/1_mapping_ont_2/minimap_t2t_sorted.merged.bam &> bam2fastq.log &
# get fasta
awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa



#母
nohup meryl count compress k=20 threads=30 memory=300 ./jin1/jin1_1.clean.fq.gz ./jin1/jin1_2.clean.fq.gz output maternal_compress.k20.meryl  > maternal_compress.log 2>&1 &
#父
nohup meryl count compress k=20 threads=20 memory=200 ./jin2/jin2_1.clean.fq.gz ./jin2/jin2_2.clean.fq.gz output paternal_compress.k20.meryl > paternal_compress.log 2>&1 &
#子
ln -s ./0_illumina/illumina/clean/*.gz .
nohup meryl count compress k=20 threads=30 memory=300 FDES220007339-1a_L1_1_clean.rd.fq.gz  FDES220007339-1a_L1_2_clean.rd.fq.gz  FDES220007339-1a_L2_1_clean.rd.fq.gz  FDES220007339-1a_L2_2_clean.rd.fq.gz output child_compress.k20.meryl > child_compress.log 2>&1 &



 meryl print ./child_compress.k20.meryl/ |head
 cd /disk202/caocy/software/merqury/
 export MERQURY=$PWD 
$MERQURY/trio/hapmers.sh maternal_compress.k20.meryl paternal_compress.k20.meryl child_compress.k20.meryl
child_compress.k30.meryl hapmers.sh maternal_compress.k30.only.meryl paternal_compress.k30.only.meryl

	 meryl print maternal_compress.k20.hapmer.meryl |head 
source /disk212/caocy/software/anaconda3/bin/activate  verkko
#/disk191_3/LabData/Pig/1_DNAseq/Third_1/0_hifi/5cell/JH_hifi.fastq.gz
#/disk191_3/Lab Data/Pig/1_DNAseq/Third_1/0_ONT3/JH-7cell/pass.fq.gz
#/disk191_3/LabData/Pig/1_DNAseq/Third_1/0_ONT/X101SC22081679-Z01-F002/20220917_0801_3E_PAM16166_3a7e2a5b/fastq_pass/onecellall.fastq.gz
nohup verkko -d ./0_verkko/ --threads 40 --local-memory 500 --local \
  --hifi JH_hifi.fastq.gz \
  --nano pass.fq.gz onecellall.fastq.gz \
  --hap-kmers ./0_verkko/meryl/maternal_compress.k20.hapmer.meryl \
              ./0_verkko/meryl/paternal_compress.k20.hapmer.meryl \
              trio > JH_verkko_3.log 2>&1 &

conda install bioconda::mbg

nohup verkko -d ./0_verkko/ --threads 40 --local-memory 500 --local \
  --hifi JH_hifi.fastq.gz \
  --nano pass.fq.gz onecellall.fastq.gz \
  --ref ./JH_t2t_6_chronly.fa \
  --hap-kmers ./0_verkko/meryl/maternal_compress.k20.hapmer.meryl \
              ./0_verkko/meryl/paternal_compress.k20.hapmer.meryl \
              trio > JH_verkko_4.log 2>&1 &


##选择更加连续的



##################挂载

###verkko
cd ./0_hifiasm_trio_UL/haplotype1/ragtag/
nohup ragtag.py scaffold -u ./JH_t2t_6_chronly.fa ./0_hifiasm_trio_UL/haplotype1/JH.hap1.fa -t 20 -o JH.haplotype1.fasta > JH.haplotype1.log 2>&1 &

cd ./0_hifiasm_trio_UL/haplotype2/ragtag
#cd ./0_hifiasm_trio_UL/haplotype2/ragtag_1/
nohup ragtag.py scaffold -u ./JH_t2t_6_chronly.fa ./0_hifiasm_trio_UL/haplotype2/JH.hap2.fa -t 20 -o JH.haplotype2.fasta > JH.haplotype2.log 2>&1 &
nohup ragtag.py scaffold -u ./JH_t2t_6_chronly.fa ./0_hifiasm_trio_UL/haplotype2/JH.hap2.fa -t 20 -o JH.haplotype2_2.fasta > JH.haplotype2.log 2>&1 &

####拆分用XY拆开的基因组进行挂载
####haplotype1用0_verkko的做骨架
nohup ragtag.py scaffold -u ./0_hifiasm_trio_UL/haplotype1/JH.haplotype1.fa ./0_verkko/assembly.haplotype1.fasta -t 20 -o JH.haplotype1.fasta > JH.haplotype1.log 2>&1 &

####haplotype1用0_hifiasm_trio_UL的做骨架
cd ./0_hifiasm_trio_UL/haplotype1/ragtag_1/
nohup ragtag.py scaffold -u ./0_hifiasm_trio_UL/haplotype_ref/JH.haplotype1.fa ./0_hifiasm_trio_UL/haplotype1/JH.hap1.fa -t 20 -o JH.haplotype1.fasta > JH.haplotype1.log 2>&1 &


####haplotype2用0_hifiasm_trio_UL的做骨架
nohup ragtag.py scaffold -u ./0_hifiasm_trio_UL/haplotype_ref/JH.haplotype2.fa ./0_hifiasm_trio_UL/haplotype2/JH.hap2.fa -t 20 -o JH.haplotype2.fasta > JH.haplotype2.log 2>&1 &

####haplotype2用0_verkko的做骨架
nohup ragtag.py scaffold -u ./0_hifiasm_trio_UL/haplotype_ref/JH.haplotype2.fa ./0_verkko/assembly.haplotype2.fasta -t 20 -o JH.haplotype2.fasta > JH.haplotype2.log 2>&1 &



cd ./0_hifiasm_trio_UL/haplotype1/JH.haplotype1.fasta/
bioawk -c fastx '{print ">"$name;print length($seq)}' ragtag.scaffold.fasta > JH_haplotype1.ragtag.summary
bioawk -c fastx '{print ">"$name;print length($seq)}' ragtag.scaffold.fasta > JH_haplotype2.ragtag.summary


#conda install -c bioconda quast
cd ./0_hifiasm_trio_UL/haplotype1/JH.haplotype1.fasta/
quast.py ./0_hifiasm_trio_UL/haplotype1/JH.haplotype1.fasta/ragtag.scaffold.fasta -o quast_results_haplotype1 &
quast.py ./0_verkko/assembly.haplotype2.fasta -o quast_results_verkko_haplotype2 &


#####分离多余的
 perl ./2_genomeAnno/split.pl ./0_hifiasm_trio_UL/haplotype1/ragtag/JH.haplotype1.fasta/ragtag.scaffold.fasta
 cat chr* > JH.haplotype1.fa
 rm chr*
 cat h1tg* > res_haplotype1.fa
 rm h1tg*
 
 perl ./2_genomeAnno/split.pl ./0_hifiasm_trio_UL/haplotype2/ragtag/JH.haplotype2.fasta/ragtag.scaffold.fasta
 cat chr* > JH.haplotype2.fa
 rm chr*
 cat h2tg* > res_haplotype2.fa
 rm h2tg*
 

###用haplotype1verrko结果用hifiasm结果补充
cd ./0_hifiasm_trio_UL/haplotype1/2_gapfill/
#--min_match 100000 --min_idy 0.8  保存paf文件
nohup tgsgapcloser --scaff ./0_hifiasm_trio_UL/haplotype1/JH.haplotype1.fasta/ragtag.scaffold.fasta --reads ./0_hifiasm_trio_UL/haplotype1/JH.hap1.fa --min_match 100000 --min_idy 0.8 --output haplotype1_1 --tgstype pb --ne --thread 36 >haplotype1_1.log 2>&1 &

###用haplotype1结果hifiasm用verrko结果补充
nohup tgsgapcloser --scaff ./0_hifiasm_trio_UL/haplotype1/ragtag_1/JH.haplotype1.fasta/ragtag.scaffold.fasta --reads ./0_verkko/assembly.haplotype1.fasta --min_match 100000 --min_idy 0.8 --output haplotype1_1 --tgstype pb --ne --thread 36 >haplotype1_1.log 2>&1 &


###用haplotype2结果hifiasm用verrko结果补充
cd ./0_hifiasm_trio_UL/haplotype2/2_gapfill/
nohup tgsgapcloser --scaff ./0_hifiasm_trio_UL/haplotype2/JH.haplotype2.fasta/ragtag.scaffold.fasta --reads ./0_verkko/assembly.haplotype2.fasta --min_match 10000 --min_idy 0.8 --output haplotype2_1 --tgstype pb --ne --thread 36 >haplotype2_1.log 2>&1 &

#verrko结果用hifiasm结果补充
nohup tgsgapcloser --scaff ./0_hifiasm_trio_UL/haplotype2/ragtag_1/JH.haplotype2.fasta/ragtag.scaffold.fasta --reads ./0_hifiasm_trio_UL/haplotype2/JH.hap2.fa --min_match 10000 --min_idy 0.8 --output haplotype2_1 --tgstype pb --ne --thread 36 >haplotype2_1.log 2>&1 &

bioawk -c fastx '{print ">"$name;print length($seq)}' haplotype2_1.scaff_seqs > JH_haplotype2.summary



##选hifiasm的结果
./0_hifiasm_trio_UL/haplotype1/3_gapfill/haplotype1_1.scaff_seqs

./0_hifiasm_trio_UL/haplotype2/2_gapfill/haplotype2_1.scaff_seqs
perl ./2_genomeAnno/split.pl ./0_hifiasm_trio_UL/haplotype1/2_gapfill/haplotype1_1.scaff_seqs
rm h1tg*
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
do
mv chr$j''_RagTag.fa chr$j.fa
#sed -i "s#chr$j''_RagTag#$j#g" chr$j.fa
cat chr$j.fa  >>haplotype1.fa
done
rm chr*
source /disk213/xieqq/miniconda/bin/activate /disk213/xieqq/miniconda/envs/BUSCO
#conda create -n BUSCO -c conda-forge -c bioconda busco=5.4.3
#source /disk212/caocy/software/anaconda3/bin/activate BUSCO

nohup busco -i ./0_hifiasm_trio_UL/haplotype1/3_gapfill/haplotype1_1.scaff_seqs -c 20 -o haplotype1 -m geno -l ./1_ASS/0_basic/0_BUSCO/mammalia_odb10/ &

nohup busco -i ./0_hifiasm_trio_UL/haplotype2/2_gapfill/haplotype2_1.scaff_seqs -c 20 -o haplotype2 -m geno -l ./1_ASS/0_basic/0_BUSCO/mammalia_odb10/  >haplotype2.log 2>&1 &

echo "
busco -i ./0_hifiasm_trio_UL/2_mapping/haplotype1.fa -c 20 -o haplotype1_chr -m geno -l ./1_ASS/0_basic/0_BUSCO/mammalia_odb10/ 
busco -i ./0_hifiasm_trio_UL/2_mapping/haplotype2.fa -c 20 -o haplotype1_chr -m geno -l ./1_ASS/0_basic/0_BUSCO/mammalia_odb10/  
busco -i ./0_hifiasm_trio_UL/haplotype2/2_gapfill/haplotype2_1.scaff_seqs -c 20 -o haplotype2 -m geno -l ./1_ASS/0_basic/0_BUSCO/mammalia_odb10/
" >> busco.sh

nohup bash busco.sh &> busco2.log &




python3 scripts/generate_plot.py –wd BUSCO_summaries
python3 generate_plot.py –wd /disk222/caocy/1_T2T/7_ASS/5_busco


#MERQURY计算QV值
#trio
source /disk212/caocy/software/anaconda3/bin/activate  python3.7
cd ./0_hifiasm_trio_UL/4_QV/
#. I have two assemblies (diploid)  Note there is no need to run merqury per-assemblies. Give two fasta files, Merqury generates stats for each and combined.
$MERQURY/merqury.sh /disk212/caocy/5_T2T/3_verkko/child.meryl /disk212/caocy/5_T2T/3_verkko/maternal.meryl /disk212/caocy/5_T2T/3_verkko/paternal.meryl ./0_hifiasm_trio_UL/2_mapping/haplotype1.fa ./0_hifiasm_trio_UL/2_mapping/haplotype2.fa out_prefix_chr &

source /disk212/caocy/software/anaconda3/bin/activate /disk212/caocy/software/anaconda3/envs/python3.7

###map
cd ./0_hifiasm_trio_UL/haplotype1/2_gapfill_3/

haplotype=/disk202/caocy/5_JH_T2T/0_hifiasm_trio_UL/2_mapping/haplotype1.fa
haplotype=/disk202/caocy/5_JH_T2T/0_hifiasm_trio_UL/2_mapping/haplotype2.fa

hifi_reads=/disk191_3/LabData/Pig/1_DNAseq/Third_1/0_hifi/5cell/JH_hifi.fastq.gz
echo "
minimap2 -ax map-pb ${haplotype} ${hifi_reads} -t 30 > aln_hifi.sam
samtools view -bS aln_hifi.sam > aln_hifi.bam
samtools sort -@ 30 -o minimap_hifi_sorted.merged.bam aln_hifi.bam
samtools index minimap_hifi_sorted.merged.bam
" >> minimap_t2t_hifi.sh

nohup bash minimap_t2t_hifi.sh >haplotype1_hifi.log 2>&1 &
nohup bash minimap_t2t_hifi.sh >haplotype2_hifi.log 2>&1 &

#samtools sort -@ 30 -o minimap_t2t_sorted.merged.bam aln_t2t.bam
#samtools index minimap_t2t_sorted.merged.bam
#ont t2t 222 samtools

ont_reads=/disk191_3/LabData/Pig/1_DNAseq/Third_1/0_ONT3/JH-7cell/pass.fq.gz 
#/disk191_3/LabData/Pig/1_DNAseq/Third_1/0_ONT3/JH-7cell/pass.fq.gz
#/disk191_3/LabData/Pig/1_DNAseq/Third_1/0_ONT3/JH-100K-7cell/pass.fq.gz 
echo "
minimap2 -ax map-ont ${haplotype} ${ont_reads} -t 36 > aln_t2t_ont.sam
samtools view -bS aln_t2t_ont.sam > aln_t2t_ont.bam
samtools sort -@ 30 -o minimap_ont_sorted.merged.bam aln_t2t_ont.bam
samtools index minimap_ont_sorted.merged.bam" >> minimap_t2t_ont.sh

nohup bash minimap_t2t_ont.sh &> minimap_t2t_ont.log &

samtools sort -@ 30 -o minimap_t2t_sorted.merged.bam aln_t2t.bam &
samtools index minimap_t2t_sorted.merged.bam &

samtools flagstat minimap.merged.bam > minimap.merged.bam.flagstat
samtools depth -aa minimap_t2t_sorted.merged.bam > depth.info


cd ./0_hifiasm_trio_UL/2_mapping/
bwa index -p haplotype1 haplotype1.fa &
bwa index -p haplotype2 haplotype2.fa &
#bioawk -c fastx '{print ">"$name;print length($seq)}' haplotype1.fa > JH_haplotype1.summary
#bioawk -c fastx '{print ">"$name;print length($seq)}' haplotype2.fa > JH_haplotype2.summary
samtools faidx haplotype1.fa
samtools faidx haplotype2.fa

haplotype=./2_mapping/haplotype1
haplotype=./2_mapping/haplotype2

#子 两个库
for i in 1 2
do
prefix=FDES220007339-1a_L${i}
echo "
bwa mem -t 30 ${haplotype} ./0_illumina/illumina/clean/${prefix}_1_clean.rd.fq.gz ./0_illumina/illumina/clean/${prefix}_2_clean.rd.fq.gz > ${prefix}.sam
samtools view -b ${prefix}.sam > ${prefix}.bam
samtools sort -@ 30 -o ${prefix}.sorted.bam ${prefix}.bam
samtools index ${prefix}.sorted.bam
rm ${prefix}.sam ${prefix}.bam
" >> map_jin3.sh
done
nohup bash map_jin3.sh &> map_jin3.log &
#############################没做
echo "
samtools merge -@ 30 jin3.bam  FDES220007339-1a_L1.sorted.bam FDES220007339-1a_L2.sorted.bam
 samtools sort -O bam -@ 30 -o jin3.sorted.bam jin3.bam
 samtools index -@ 30 jin3.sorted.bam
" >> map_jin3_M.sh

nohup bash map_jin3_M.sh &> map_jin3_M.log &
 
samtools flagstat jin3.sorted.bam > jin3.sorted.bam.flagstat
 
 
 
bedtools makewindows -g ./0_hifiasm_trio_UL/2_mapping/haplotype1.fa.fai -w 5000000 >5mb.bed
#bedtools nuc -fi ./0_hifiasm_trio_UL/2_mapping/haplotype1.fa -bed 50mb.bed | cut -f 1-3,5 > 5mb.bed &

nohup mosdepth -b 5mb.bed hifi_depth ./0_hifiasm_trio_UL/2_mapping/haplotype1/minimap_hifi_sorted.merged.bam &
nohup mosdepth -b 5mb.bed ont_depth ./0_hifiasm_trio_UL/2_mapping/haplotype1/minimap_ont_sorted.merged.bam &

#nohup mosdepth -b 5mb.bed ngs_depth ./1_ASS/1_mapping_ngs/2_bam/jin3.sorted.bam &

bedtools makewindows -g ./0_hifiasm_trio_UL/2_mapping/haplotype2.fa.fai -w 5000000 >5mb.bed
#bedtools nuc -fi ./0_hifiasm_trio_UL/2_mapping/haplotype1.fa -bed 50mb.bed | cut -f 1-3,5 > 5mb.bed &

nohup mosdepth -b 5mb.bed hifi_depth ./0_hifiasm_trio_UL/2_mapping/haplotype2/minimap_hifi_sorted.merged.bam &
nohup mosdepth -b 5mb.bed ont_depth ./0_hifiasm_trio_UL/2_mapping/haplotype2/minimap_ont_sorted.merged.bam &

 
 

 
#################################################################################
##两个基因组间的snp  indel
cd ./0_hifiasm_trio_UL/3_comparing/
echo "
minimap2 -cx asm5 -t8 --cs ./0_hifiasm_trio_UL/2_mapping/haplotype1.fa ./0_hifiasm_trio_UL/2_mapping/haplotype2.fa > asm.paf   &
sort -k6,6 -k8,8n asm.paf > asm.srt.paf           # sort by reference start coordinate
k8 paftools.js call asm.srt.paf -f ./0_hifiasm_trio_UL/2_mapping/haplotype1.fa -s D > out.vcf
paftools.js call asm.srt.paf -f ./0_hifiasm_trio_UL/2_mapping/haplotype1.fa -s D -L 20000 > out.vcf" >> minimap_GvG.sh

nohup bash minimap_GvG.sh &> minimap_GvG.log &


