# Structural Variants


minimap2 -a -x asm5 --cs -r2k -t 10 ./index/Sscrofa11.fa ./JH_t2t_all.fa > alignments_t2t.sam \  
samtools sort -m4G -@4 -o  alignments_t2t.sorted.bam alignments_t2t.sam \     # sort by reference start coordinate
samtools index alignments_t2t.sorted.bam \
svim-asm haploid ./8_SV alignments_t2t.sorted.bam ./index/Sscrofa11.fa

./envs/svimasm_env/bin/svim-asm haploid /./S_T alignments_S_T.sorted.bam ./Sscrofa11_1.fa

./ensembl-vep/vep -i variants.vcf  -o t2t_sv.vep --species sus_scrofa  --cache --dir_cache ./.vep/  --sift b  --symbol  --variant_class


svim-asm haploid ./1_T2T/8_SV/3_reft2t/ ./1_T2T/8_SV/0_t2t/alignments_t2t.sorted.bam /disk222/caocy/1_T2T/7_ASS_2/0_basic/t2t_fa_all/JH_t2t_all.fa

minimap2 -a -x asm5 --cs -r2k -t 10  ./1_T2T/7_ASS_2/0_basic/t2t_fa_all/JH_t2t_all.fa ./1_T2T/7_ASS/1_mapping_ngs/index/Sscrofa11.fa > alignments_t2t.sam   
samtools sort -m4G -@4 -o  alignments_t2t.sorted.bam alignments_t2t.sam    # sort by reference start coordinate
samtools index alignments_t2t.sorted.bam 
svim-asm haploid ./1_T2T/8_SV/3_reft2t/ alignments_t2t.sorted.bam ./1_T2T/7_ASS_2/0_basic/t2t_fa_all/JH_t2t_all.fa

svim-asm haploid ./1_T2T/8_SV/0_t2t/50_250bp/ alignments_t2t.sorted.bam ./1_T2T/7_ASS_2/0_basic/t2t_fa_all/JH_t2t_all.fa \
--min_sv_size 50 --max_sv_size 250

