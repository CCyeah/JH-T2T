##Genome size estimation

nohup jellyfish count -m 17 -s 10G -t 12 -C -o kmercount JH.fastq &

jellyfish stats kmercount -o kmer_count
jellyfish histo -o K20.histo K20

kmercount -o 17merFreq
kmercount -o jelly.log

### convert bam to fastq

nohup bam2fastq -o JH_hifi ./JH/JHpart1/hifi/m64270e_220824_022211.hifi_reads.bam ./JH/JHpart1/hifi/m64284e_220824_022534.hifi_reads.bam ./JH/JHt2t/Data-X101SC22081679-Z01-J005/ZJUHiC-1/Sequel-IIe/FPAC220011846-1A/m64168e_220917_191220.hifi_reads.bam  ./JH/JHt2t/Data-X101SC22081679-Z01-J005/ZJUHiC-1/Sequel-IIe/FPAC220011846-1A-1/m64168e_220919_012408.hifi_reads.bam &> bam2fastq.log &

#5个cell合并成一个文件
cd /disk195/caocy/5_genome/T2T/0_hifi
nohup bam2fastq -o JH_hifi ./m64270e_220928_013944.hifi_reads.bam \
./m64270e_220824_022211.hifi_reads.bam  ./m64284e_220824_022534.hifi_reads.bam \
./m64168e_220917_191220.hifi_reads.bam  ./m64168e_220919_012408.hifi_reads.bam &> bam2fastq.log &

### 5cell
nohup hifiasm -o JH_hifi.asm -t32 ./JH_hifi.fastq.gz 2> JX_hifi.log &
awk '/^S/{print ">"$2;print $3}' JH_hifi.asm.bp.p_ctg.gfa > JH_hifi.fa &


#ont
##QC
nohup cat ./80k.onecellall.fastq.gz ./jinhua/pass.fq.gz   ./JH-100K-PAM75762/pass.fq.gz ./JH-100K-PAM77283/pass.fq.gz ./JH-100K-PAM75878/pass.fq.gz  >merge_all_final.fastq.gz &

nohup ./seq_stat -g 2660.49M -d 6 merge.fofn > merge.log 2>&1 &
#gunzip 管道处理压缩数据
nohup gunzip -c ./fastq_pass/onecellall.fastq.gz | NanoFilt -q 7 -l 80000  | gzip > onecellall.NanoFilt.fastq.gz  &> NanoFilt.log &

nohup NanoPlot --fastq ./pass.fq.gz  -o fastq-plots  -t 10  --plots hex dot &

nohup nextDenovo  run.cfg  &


###hifi
cd /disk222/caocy/1_T2T/5_gapfill_TGS_2/1_hifi/
nohup tgsgapcloser --scaff /disk222/caocy/1_T2T/4_gapfilling_M/4_gapfilling_4/genome.review4.assembly.FINAL.fasta \
--reads /disk222/caocy/1_T2T/5_gapfill_TGS_2/JH_hifi.fasta \
        --output hifi_fill_outfile \
        --thread 48 \
        --min_idy 0.5 --min_match 2000 --ne \
        >pipe.log 2>pipe.err &

###################################################################################ONT
nohup tgsgapcloser --scaff /disk191_3/LabData/Pig/1_DNAseq/Third_1/1_ASS/X101SC22081679-Z01-F005/genome.review4.assembly.FINAL.fasta \
      --reads  /disk222/caocy/1_T2T/5_gapfill_TGS_2/nd.asm.fasta \
      --min_match 2000 \
	  --min_idy 0.6 \
	  --thread 36 \
      --output result_6  > ont_all_fill.log 2>ont_all_fill.err &







