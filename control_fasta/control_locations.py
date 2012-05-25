#### test on older data..
#### test new info...
from random_noncoding_region import recursive_merge_both, non_coding
from Bio.Seq import Seq
from itertools import chain
from collections import deque, defaultdict
from flatfeature import Bed
from pyfasta import Fasta
import csv

def get_genespace(cns_file):
    """given a cns list gets genespace for assigned gene"""
    d = defaultdict(list)
    for line in open(cns_file):
        if line[0] == '#': continue
        #cns_id,seqid,start,end,strand,accn=line.strip().split('\t')[:6]
        cns_id,qaccn,seqid,start,end,strand,saccn,sseqid,sstart,ssend,sstrand,ev,url=line.strip().split(',')
        
        start,end,sstart,ssend = map(int,[start,end,sstart,ssend])
        if start > end:
            start,end = end,start
        #if sstart > ssend:
        #    ssend,sstart = sstart,ssend
        d[qaccn].append((start,end))
        #d[saccn].append((sstart,ssend))
    return d

def get_pos(gene,locs,cnss):
    postion_info = {}
    ########### postion info #############
    genespace_start, genespace_end = cnss[0][0],cnss[-1][1]
    postion_info['3_utr'] = (locs[-1][1],int(gene['end']))
    postion_info['5_utr'] = (int(gene['start']),locs[0][0])
    postion_info['intronic'] = (locs[0][0],locs[-1][1])
    postion_info['3_prox'] = (int(gene['end']) ,int(gene['end'])+ 3000 )
    postion_info['5_prox'] = (int(gene['start']) - 3000,int(gene['start']))
    postion_info['5_distal'] = (min(genespace_start, int(gene['start']) - 3000),int(gene['start']) - 3000)
    postion_info['3_distal'] =(int(gene['end']) + 3000 ,max(int(gene['end']) +3000, genespace_end))
    if gene['strand'] == '-':
        postion_info['3_utr'],postion_info['5_utr'],postion_info['3_prox'],postion_info['5_prox'],postion_info['5_distal'],postion_info['3_distal'] = postion_info['5_utr'],postion_info['3_utr'],postion_info['5_prox'],postion_info['3_prox'],postion_info['3_distal'],postion_info['5_distal']

    ###### something like this double check for strand info ######
    return postion_info


def write_to_pos_fasta(bed,accn,locs,cnsspace,fhs,f):
    postion_info = get_pos(accn,locs,cnsspace)
    remove_cds(bed,accn,cnsspace[0][0],cnsspace[-1][1],postion_info,fhs,f)

def open_files(postions):
    fhs = {}
    for postion in postions:
        fh = open('{0}_control.fasta'.format(postion),'wb')
        fhs[postion] = fh
    return fhs

def main(cns_file,bedpath,fastapath):
    genespace = get_genespace(cns_file)
    bed = Bed(bedpath)
    f = Fasta(fastapath)
    handles = ['3_utr','5_utr','intronic','5_prox','5_distal','3_prox','3_distal']
    fhs = open_files(handles)
    for gene in genespace.keys():
        #cnsspace = genespace[gene]
        try:
            accn = bed.accn(gene)
        except KeyError: continue
        cnsspace = [(max(0,accn['start'] - 12000), accn['end'] + 12000)]
        #print "GENESPACE {0}".format(cnsspace)
        locs = accn['locs']
        locs.sort()
        cnsspace.sort()
        write_to_pos_fasta(bed,accn,locs,cnsspace,fhs,f)
    
def remove_cds(bed,gene,search_start,search_end,postion_info,fhs,f):
    feats = bed.get_features_in_region(gene["seqid"],search_start,search_end)
    #print "feats ...{0}".format(feats)
    locs =[feat["locs"] for feat in feats]
    #print search_start,search_end
    ### looks for cns masked regions and cds
    flatten_locs = chain.from_iterable(locs)
    flatten_locs = list(flatten_locs)
    ##iistart_stops =[(int(start),int(end)) for (start,end) in locs]
    start_stops =[(int(start),int(end)) for (start,end) in flatten_locs]
    start_stops_in_range = []
    for s_start,s_stop in start_stops:
        if s_start <= search_end and s_stop >= search_start:
            start_stops_in_range.append((s_start,s_stop))
    start_stops_in_range.sort()
    if len(start_stops_in_range) > 0:
        if len(start_stops_in_range) > 1:
            start_stops_m = recursive_merge_both(deque(start_stops_in_range))
        else:
            start_stops_m = start_stops_in_range
        merged_feats = list(start_stops_m)
        merged_feats.sort()
        off_by_one = non_coding(merged_feats,search_start,search_end)
        #print off_by_one
    else:
        off_by_one = [(search_start,search_end)]
    get_fasta(gene,off_by_one,f,postion_info,fhs)


def get_fasta(accn,off_by_one,f,pos_info,fhs):
    off_by_one.sort()
    seq = f[accn['seqid']]
    n = 0
    strand = accn['strand']
    for start,stop in off_by_one:
        fasta = seq[start-1:stop-1]
        if len(fasta) <  15:
            print start,stop,accn['accn']
            continue
        n += 1
        #print seq[start-1:stop-1]
        pos = [pos for pos in pos_info.keys() if start in range(pos_info[pos][0],pos_info[pos][1])]
        #print pos,start,stop,pos_info.values()
        new_fasta = fhs[pos[0]]
        w = ">cns{1}_{0}\n".format(accn['accn'],n)
        seq_w = "{0}\n".format(fasta)
        new_fasta.write(w)
        new_fasta.write(seq_w)
    
#main('/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_golden/golden_cns_from_dropbox.csv','/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_golden/thaliana_v8.with_new_cns_mask.bed','/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/thaliana_v8.fasta')
#main('/Users/gt/Desktop/paper/G-box-seq/rice_rice/rice_j_rice_j.cns.assigned.csv.local','/Users/gt/Desktop//paper/G-box-seq/rice_rice/rice.with_new_cns_mask.bed','/Users/gt/Desktop/paper/G-box-seq/rice_rice/rice_j.fasta')
#main('/Users/gt/Desktop/paper/G-box-seq/set_set/setaria_n_setaria_n.cns.assigned.csv.local','/Users/gt/Desktop/paper/G-box-seq/set_set/set.with_new_cns_mask.bed','/Users/gt/Desktop/paper/G-box-seq/set_set/setaria_n.fasta')
#main('/Users/gt/Desktop/paper/G-box-seq/rice_set/rice_j_setaria_n.cns.assigned_real.csv.local','/Users/gt/Desktop/paper/G-box-seq/rice_set/rice_j_set.nolocaldups.with_new_cns_mask.bed','/Users/gt/new/rice_j.fasta')
main('/Users/gt/Desktop/paper/G-box-seq/rice_sorg/rice_j_sorghum_nn.cns.assigned_real.csv.local','/Users/gt/Desktop/paper/G-box-seq/rice_sorg/rice_j_sorg.nolocaldups.with_new_cns_mask.bed','/Users/gt/new/rice_j.fasta')

