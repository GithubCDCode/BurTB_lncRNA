#get clean
fastp --qualified_quality_phred 5 --unqualified_percent_limit 50  --n_base_limit 15 --min_trim_length 10 --overlap_len_require 30 --overlap_diff_limit 1 --overlap_diff_percent_limit 10  --length_required 150 --length_limit 150 --trim_poly_g --thread 1 --html qc_report.html -j sample.json -i sample_1.fq.gz -I sample_2.fq.gz -o sample_1.clean.fq.gz -O sample_2.clean.fq.gz
#mapping
hisat2 -x genome_index -p 4 --dta -t --phred33 --rna-strandness RF -1 sample_1.clean.fq.gz -2 sample_2.clean.fq.gz --un-conc-gz sample.unmap.fq.gz 2> align.log | samtools sort -O BAM --threads 4 -o sample.bam -
#assemble
stringtie sample.bam -p 4 -G genome.gtf -o sample.gtf --rf
stringtie --merge -G genome.gtf -l novel -o merge.gtf sample1.gtf sample2.gtf sample3.gtf <...>
gffcompare -R -r genome.gtf -o ./gffcompare merge.gtf
#filter
cuffquant --library-type fr-firststrand -p 2 --max-bundle-frags 1165754 -o ./sample/sample.cuffquant novel.merge.gtf sample.bam
cuffnorm --library-type fr-firststrand --num-threads 2 -o ./cuffnorm -L sample1,sample2,<...> novel.merge.gtf ./sample1/abundances.cxb ./sample2/abundances.cxb <...>
根具isoforms.fpkm_table中表达量筛选
pfam_scan.pl -fasta protein.fa<表达量筛选得到的序列> -dir database/pfam -outfile pfam_scan.out -cpu 4
CNIT.py -f novel.exp_filter.fasta -m ve -o ./
CPC2.py -i novel.exp_filter.fasta -o CPC2.result.xls
pfam CNIT CPC2结果取交集
#quant
stringtie -e -p 8 -G all_transcripts.gtf -A genes.xls -o transcript.gtf sample.bam
#diff 
edgeR --count genes.readcount.xls -condition condition.xls --cpname group1vsgroup2 --fc 1 --padj 0.05 --outdir ./

