import argparse
import sys
from os import path
from os import system
from os import chdir
from HTSeq import GFF_Reader
from collections import defaultdict
import json

parser=argparse.ArgumentParser(description = 'generate co_location table')
parser.add_argument('--mRNA_gtf', help = 'mRNA gtf file.', required = True)
parser.add_argument('--lncRNA_gtf', help = 'lncRNA or TUCP gtf file .', required = True)
parser.add_argument('--name', help = 'transcript type.',choices = ['lncRNA','TUCP'], required = True)
parser.add_argument('--distance', help = 'co_location distance.', default = "100000")
parser.add_argument('--tr_info', help = 'JSON format transcript info file.')
parser.add_argument('--gtf2bed', help = 'gtf2bed',required = True)
parser.add_argument('--bedtools', help = 'bedtools',required = True)
parser.add_argument('--out_dir', help = 'output dir',required = True)
argv = vars(parser.parse_args())

def generate_dir(dir) :
    if not path.exists(dir) :
        assert not system("mkdir {dir}".format(**locals()))
    else :
        pass
    return 0

def get_info_from_gtf(gtf) :
    gtf_info = {}
    for eachline in GFF_Reader(gtf):
        gene_id   = eachline.attr['gene_id']
        if 'transcript_type' in eachline.attr :
            gene_type = eachline.attr['transcript_type']
        elif 'gene_biotype' in eachline.attr :
            gene_type = eachline.attr['gene_biotype']
        elif 'gene_type' in eachline.attr :
            gene_type = eachline.attr['gene_type']
        else :
            print ("No gene type information in gtf.!")
            sys.exit(0)
        genename = '--'
        if 'gene_name' in eachline.attr :
            genename = eachline.attr['gene_name']
        chrom  = eachline.iv.chrom
        start  = eachline.iv.start
        end    = eachline.iv.end
        strand = eachline.iv.strand
        if gene_id in gtf_info :
            if start < gtf_info[gene_id]['start'] :
                gtf_info[gene_id]['start'] = start
            if end > gtf_info[gene_id]['end'] :
                gtf_info[gene_id]['end'] = end
        else :
            gtf_info.setdefault(gene_id,{})["gene_id"] = eachline.attr['gene_id']
            gtf_info.setdefault(gene_id,{})["chrom"] = chrom
            gtf_info.setdefault(gene_id,{})["start"] = start
            gtf_info.setdefault(gene_id,{})["end"] = end
            gtf_info.setdefault(gene_id,{})["strand"] = strand
            gtf_info.setdefault(gene_id,{})["tr_type"] = gene_type
            gtf_info.setdefault(gene_id,{})["gene_name"] = genename
    return gtf_info

def get_tr_gene_dict(gtf):
    tr_gene_dict = {}
    for eachline in GFF_Reader(gtf):
        gene_id = eachline.attr['gene_id']
        transcript_id = eachline.attr['transcript_id']
        if transcript_id not in tr_gene_dict :
            tr_gene_dict[transcript_id] = gene_id
    return tr_gene_dict

def get_up_down_info(target_strand,target_start,target_end,query_start,query_end) :
    distance = 0
    up_down_stat = None
    if target_start < query_start :
        distance = query_start - target_end
        if target_strand == "+" :
            up_down_stat = "downstream"
        else :
            up_down_stat = "upstream"
    else :
        distance = target_start - query_end
        if target_strand == "+" :
            up_down_stat = "upstream"
        else :
            up_down_stat = "downstream"
    return distance,up_down_stat

def generate_info(gtf_dict,query_id) :
    query_start = gtf_dict[query_id]['start']
    query_end = gtf_dict[query_id]['end']
    query_strand = gtf_dict[query_id]['strand']
    query_info = "{query_start}\t{query_end}\t{query_strand}".format(**locals())
    return query_info

#gtf2bed = "/TJPROJ1/RNA_R/pipline/ncRNA/V3.2.0_1211/V3.2.2/bin/gtf2Bed.pl"
#cuffcompare = "/PUBLIC/software/RNA/cufflinks-2.1.1/cuffcompare"
#bedtools = "/PUBLIC/software/HUMAN/bin/bedtools"


mRNA_gtf = argv['mRNA_gtf'].strip()
lncRNA_gtf = argv['lncRNA_gtf'].strip()
distance = argv['distance'].strip()
out_dir = argv['out_dir'].strip()
name = argv['name'].strip()
tr_info = argv['tr_info'].strip()
gtf2bed = argv['gtf2bed'].strip()
bedtools = argv['bedtools'].strip()
generate_dir(out_dir)
gene_info_dict=get_info_from_gtf(mRNA_gtf)
lncRNA_info_dict=get_info_from_gtf(lncRNA_gtf)
if __name__ == '__main__':
    # check file and directory
    for each in [mRNA_gtf,lncRNA_gtf,out_dir] :
        if not path.exists(each) :
            print ("{each} does no exist !".format(**locals()))
            sys.exit(0)

#    bedtools get mRNA overlapping ,upstream and downstream lncRNA
    if not path.exists('{out_dir}/mRNA.bed'.format(**locals())) :
        system('perl {gtf2bed} {mRNA_gtf} gene > {out_dir}/mRNA.bed'.format(**locals()))
    system('perl {gtf2bed} {lncRNA_gtf} transcript > {out_dir}/{name}.transcript.bed'.format(**locals()))
    system('{bedtools} intersect -a {out_dir}/mRNA.bed -b {out_dir}/{name}.transcript.bed -wo > {out_dir}/mRNA_overlapping_{name}.txt'.format(**locals()))    
    system('{bedtools} window -a {out_dir}/mRNA.bed -b {out_dir}/{name}.transcript.bed -l {distance} -r {distance} > {out_dir}/mRNA_up_down_stream_{name}.txt'.format(**locals()))

    ## get mRNA gene and lncRNA transcript information
    lncRNA_mRNA_pairs_dict = {}
    #lncRNA_mRNA_gene_pairs_dict = {}
    with open(tr_info,'r') as tr_info_file :
        tr_info_dict = json.load(tr_info_file)
    ## overlapping information
    overlapping_info_file = path.join(out_dir,'mRNA_overlapping_{name}.txt'.format(**locals()))
    with open(overlapping_info_file,'r') as overlapping_info :
        overlapping_info = [line.strip() for line in overlapping_info]
        for eachline in overlapping_info :
            each_overlapping_info = eachline.split('\t')
            mRNA_id = each_overlapping_info[3]
            #mRNA_gene = tr_info_dict[mRNA_id]['gene_id']
            mRNA_strand = each_overlapping_info[5]
            lncRNA_id = each_overlapping_info[15]
            lncRNA_strand = each_overlapping_info[17]
            overlapping_stat = '(%s)' % each_overlapping_info[24]
            if mRNA_strand == lncRNA_strand :
                location = 'sense'
            else :
                location = 'antisense'
            lncRNA_mRNA_pairs_dict.setdefault(lncRNA_id,{})[mRNA_id] = [overlapping_stat,location]
            #lncRNA_mRNA_gene_pairs_dict.setdefault(lncRNA_id,{})[mRNA_gene] = 1
    ## UP Down information
    up_down_stream_lncRNA_file = '{out_dir}/mRNA_up_down_stream_{name}.txt'.format(**locals())
    with open(up_down_stream_lncRNA_file,"r") as up_down_info :
        for eachline in up_down_info :
            if eachline.strip() != "" :
                each_up_down_info = eachline.strip().split("\t")
                lncRNA_id = each_up_down_info[15]
                mRNA_id = each_up_down_info[3]
                #mRNA_gene = tr_info_dict[mRNA_id]['gene_id']
                mRNA_strand = each_up_down_info[5]
                mRNA_start = int(each_up_down_info[6])
                mRNA_end = int(each_up_down_info[7])                
                lnRNA_start = int(each_up_down_info[18])
                lncRNA_end = int(each_up_down_info[19])
                if lncRNA_id in lncRNA_mRNA_pairs_dict and mRNA_id in lncRNA_mRNA_pairs_dict[lncRNA_id] :
                    continue
                distance,up_down_stat = get_up_down_info(mRNA_strand,mRNA_start,mRNA_end,lnRNA_start,lncRNA_end)
                lncRNA_mRNA_pairs_dict.setdefault(lncRNA_id,{})[mRNA_id] = [distance,up_down_stat]
                #lncRNA_mRNA_gene_pairs_dict.setdefault(lncRNA_id,{})[mRNA_gene] = 1
    ## output lncRNA or TUCP co_located mRNA

    lncRNA_co_location_file = "{out_dir}/{name}_mRNA_co_location_results.xls".format(**locals())
    lncRNA_pairs_file = "{out_dir}/{name}_mRNA_co_location_pairs.json".format(**locals())
    lncRNA_co_location_results = open(lncRNA_co_location_file,"w")  
    lncRNA_co_location_results.write("Chromosome\t{name}_ID\t{name}_Gene_ID\t{name}_Gene_Symbol\t{name}_Status\t{name}_Start\t{name}_End\t{name}_Strand\tGene_Id\tGene_Symbol\tGene_Start\tGene_End\tGene_Strand\tDistance\tLocation\n".format(**locals()))
    for each_lncRNA in lncRNA_mRNA_pairs_dict :
        Chromosome = tr_info_dict[each_lncRNA]['chrom']
        lncRNA_info = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (each_lncRNA,tr_info_dict[each_lncRNA]['gene_id'],tr_info_dict[each_lncRNA]['gene_name'],tr_info_dict[each_lncRNA]['status'],tr_info_dict[each_lncRNA]['start'],tr_info_dict[each_lncRNA]['end'],tr_info_dict[each_lncRNA]['strand'])
        for each_mRNA in lncRNA_mRNA_pairs_dict[each_lncRNA] :
            mRNA_id = each_mRNA
            mRNA_info = '%s\t%s\t%s\t%s\t%s' % (mRNA_id,gene_info_dict[mRNA_id]['gene_name'],gene_info_dict[mRNA_id]['start'],gene_info_dict[mRNA_id]['end'],gene_info_dict[mRNA_id]['strand'])
            Distance = lncRNA_mRNA_pairs_dict[each_lncRNA][each_mRNA][0]
            Location = lncRNA_mRNA_pairs_dict[each_lncRNA][each_mRNA][1]
            results = "{Chromosome}\t{lncRNA_info}\t{mRNA_info}\t{Distance}\t{Location}\n".format(**locals())                                        
            lncRNA_co_location_results.write(results)
    lncRNA_co_location_results.close()
    #with open(lncRNA_pairs_file,"w") as lncRNA_pairs :
    #    json.dump(lncRNA_mRNA_gene_pairs_dict,lncRNA_pairs)


