#!/usr/bin/python
# updated by szchen

import os
import sys
import argparse
import json
from HTSeq import GFF_Reader
from scipy.stats import pearsonr 


def get_id_info(tr_id,tr_info_dict,tr_type = 'mRNA'):
    tr_info = ''
    if tr_id in tr_info_dict :
        tr_gene_id = tr_info_dict[tr_id]['gene_id']
        tr_gene_name = tr_info_dict[tr_id]['gene_name']
        tr_gene_type = tr_info_dict[tr_id]['gene_type']
        if tr_type == 'mRNA' :
            tr_info = '{tr_id}\t{tr_gene_id}\t{tr_gene_name}'.format(**locals())
        else :
            tr_status = tr_info_dict[tr_id]['status']
            tr_info = '{tr_id}\t{tr_gene_id}\t{tr_gene_name}\t{tr_status}'.format(**locals())
    else :
        print (' {tr_id} not in transcript info file !'.format(**locals()))
        sys.exit(0)
    return tr_info,tr_gene_id


def get_info_from_gtf(gtf) :
    gtf_info = {}
    for eachline in GFF_Reader(gtf):
        gene_id = eachline.attr['gene_id']
        genename = '--'
        if 'gene_name' in eachline.attr :
            genename = eachline.attr['gene_name']
        chrom = eachline.iv.chrom
        start = eachline.iv.start
        end = eachline.iv.end
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
            gtf_info.setdefault(gene_id,{})["gene_name"] = genename
    return gtf_info


def read_file(file):
    with open(file,'r') as df:
        for line in df:
            yield line


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prediction of lncRNA\'s trans function.\nR function "cor.test" is applied to calculate the cor.')
    parser.add_argument('--lnc',required=True,help='Expression table of lncRNA, FPKM or RPKM.')
    parser.add_argument('--mrna',required=True,help='Expression table of mRNA, FPKM or RPKM.')
    parser.add_argument('--mRNA_gtf', help = 'mRNA gtf file.', required = True)
    parser.add_argument('--prefix',required=True,help='Output file prefix.')
    parser.add_argument('--number',required=True,help='Output file split number.')
    parser.add_argument('--tr_info',required=True,help='Transcript_info, including chromosome, position, transcript_id,gene_id,gene_name and gene description.(json dict)')
    parser.add_argument('--outdir',required=True,help='The output dir, 3 files produced: lncRNA_mRNA.cortest.xls, lncRNA_mRNA.trans.xls, lncRNA_mRNA.trans.list.')
    parser.add_argument('--corr',default=0.95,type=float,help='The threshold of correlation coefficiency, default 0.95.')
    parser.add_argument('--relation',default='both',choices=['both','negative','positive'],help='Correlation type between lncRNA and mRNA, negative, positive or both, default both.')
    argv=parser.parse_args()

    out_prefix = argv.prefix.strip()
    output_dir = argv.outdir.strip()
    split_number = argv.number.strip()
    transcript_info_file = argv.tr_info.strip()
    
    if not os.path.exists(output_dir):
        os.system('mkdir -p '+output_dir)

    output_name = '{out_prefix}_mRNA_cortest_{split_number}.xls'.format(**locals())
    co_exp_name = '{out_prefix}_mRNA_co_expression_{split_number}.json'.format(**locals())
    co_exp_pair_name = '{out_prefix}_mRNA_pairs_{split_number}.json'.format(**locals())
    output_file = os.path.join(output_dir,output_name)
    co_exp_file = os.path.join(output_dir,co_exp_name)
    co_exp_pair_file = os.path.join(output_dir,co_exp_pair_name)
    #lnc_mrna_cor = open(output_file,'w')

    with open(transcript_info_file,"r") as transcript_info_file_info :
        transcript_info_dict = json.load(transcript_info_file_info)

    mRNA_gtf = argv.mRNA_gtf.strip()
    gene_info_dict = get_info_from_gtf(mRNA_gtf)

    co_express_dict = {}
    co_express_pair_dict = {}
    j = 1
    trans = {}

    for line in read_file(argv.lnc):
        array_lnc = line.strip().split('\t')
        lnc_rpkm = list(map(float, array_lnc[1:]))
        if sum(lnc_rpkm) < 0.4:
            continue
        lncRNA_id = array_lnc[0]
        j = 1
        for each in read_file(argv.mrna):
            if j == 1:
                j += 1
                continue
            array_mrna = each.strip().split('\t')
            pc_rpkm = list(map(float, array_mrna[1:]))
            if sum(pc_rpkm) < 0.4:
                continue
            mRNA_id = array_mrna[0]
            pe, pv = pearsonr(lnc_rpkm, pc_rpkm)
            if abs(pe) > argv.corr:                        
                lncRNA_info,lncRNA_gene = get_id_info(array_lnc[0],transcript_info_dict,out_prefix)
                mRNA_info = gene_info_dict[mRNA_id]['gene_name']
                mRNA_gene = gene_info_dict[mRNA_id]['gene_id']
                cor_stat = str(pe)
                cor_pval = str(pv)
                co_express_output = '{lncRNA_info}\t{mRNA_id}\t{mRNA_info}\t{cor_stat}\t{cor_pval}\n'.format(**locals())
                co_express_dict.setdefault(lncRNA_id,{})[mRNA_id] = co_express_output
                co_express_pair_dict.setdefault(lncRNA_id,{})[mRNA_gene] = 1

    with open(co_exp_file,'w') as co_exp_file_json :
        json.dump(co_express_dict,co_exp_file_json)

    with open(co_exp_pair_file,'w') as co_exp_pair_json :
        json.dump(co_express_pair_dict,co_exp_pair_json)
