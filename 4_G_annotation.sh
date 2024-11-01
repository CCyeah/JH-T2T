##
#liftoff -g {REF_GFF} -o {output} -p 10 -dir 02_anno/{wildcards.sample}_intermediate_files {input} {INDEX_REF}
nohup liftoff -g Sus_scrofa.Sscrofa11.1.109.chr.gff3 -o t2t.gene.gff3 -p 10 result.scaff_seqs Sus_scrofa.Sscrofa11.1.dna.toplevel.fa &
wget https://ftp.ensembl.org/pub/release-109/gff3/sus_scrofa/Sus_scrofa.Sscrofa11.1.109.chr.gff3.gz
gunzip Sus_scrofa.Sscrofa11.1.109.chr.gff3.gz

samtools faidx Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
sed -i "s#Chr19#ChrX#g" JH_t2t.fa
sed -i "s#Chr20#ChrY#g" JH_t2t.fa
samtools faidx JH_t2t.fa

nohup liftoff -g Sus_scrofa.Sscrofa11.1.109.chr.gff3\
 -o t2t.gene.gff3 -p 10 JH_t2t.fa Sus_scrofa.Sscrofa11.1.dna.toplevel.fa &


fin=./1_T2T/7_ASS_2/0_basic/t2t_fa/JH_t2t.fa
#T2T=./1_T2T/7_ASS_2/0_basic/t2t_fa/JH_t2t.fa
BuildDatabase -engine ncbi -name T2T ${fin}
#nohup RepeatModeler -database NRS -pa 10 -LTRStruct >& PanRepet.out &
nohup RepeatModeler -database T2T -pa 10 >& T2TRepet.out &
cat ./PanGenome/software/repeatDatabase/Mammalias.fasta T2T-families.fa > pigt2t_repeats.fasta
# repeatMasker 注释重复序列
#nohup RepeatMasker -xsmall -poly -pa 20 -lib ./1_T2T/6_Anno/1_repeatmodel/pigt2t_repeats.fasta -gff -engine ncbi ./1_T2T/5_gapfill_TGS/17_JLcontig/result.scaff_seqs &
# nohup RepeatMasker -xsmall -poly -pa 5 -lib Mammalias.fasta -gff -engine ncbi nonRefSeq_qc2.fasta &
nohup RepeatMasker -xsmall -poly -pa 20 -lib ./1_T2T/6_Anno_2/0_repeat_model/pigt2t_repeats.fasta -gff -engine ncbi ./1_T2T/7_ASS_2/0_basic/t2t_fa/JH_t2t.fa &

source ./software/anaconda3/bin/activate ./software/anaconda3/envs/T2T


liftofftools all -r susScr11.fa \
-t result.scaff_seqs \
-rg Sus_scrofa.Sscrofa11.1.108.gff3 \
-tg JH_t2t.gff &

#bioawk -c fastx '{print ">"$name;print length($seq)}'  susScr11.fa > ref.summary

# extract fasta sequence

gffread -w transcripts_all.fa -g ./1_T2T/6_Anno_2/JH_t2t.fa ./2_genomeAnno/1_RNA/merge.gtf
##检查


#AACCCT AGGGTT
#AGGGTT TCCCAA 
#conda install -c bioconda tidk
#tidk search --string CCCTAAA  --output try.tsv --dir ./1_T2T/6_Anno_2/2_Te/ --extension tsv ./1_T2T/7_ASS_2/0_basic/t2t_fa/JH_t2t.fa
tidk search --string CCCTAA \
 --output onerepeat.tsv --dir ./1_T2T/6_Anno_2/2_Te/ --extension tsv ./1_T2T/7_ASS_2/0_basic/t2t_fa/JH_t2t.fa

tidk search --string TTAGGG \
 --output onerepeat_R.tsv --dir ./1_T2T/6_Anno_2/2_Te/ --extension tsv ./1_T2T/7_ASS_2/0_basic/t2t_fa/JH_t2t.fa

tidk search --string CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA \
 --output onerepeat_45n.tsv --dir ./1_T2T/6_Anno_2/2_Te/ --extension tsv ./1_T2T/7_ASS_2/0_basic/t2t_fa/JH_t2t.fa

