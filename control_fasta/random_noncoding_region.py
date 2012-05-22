import random
from flatfeature import Bed
from pyfasta import Fasta
from collections import deque
from itertools import chain
from Bio.Seq import Seq
"""input gene space of intrest and features of intrest in space
output list of random features in the gene space same size as other feature not
overlapping other ones"""

# add cns to bed file
# mask genome noncodding seq...

def gene_size_dict(tab_sep_file):
    """select accn, group_concat(qdiff) from
    rice_b_setaria64_gturco_2011_4_11app_real group by accn"""
    handle = open(tab_sep_file)
    fh = handle.read()
    gene_size_dict = {}
    for line in fh.split("\n")[:-1]:
        gene, cnss = line.split("\t")
        gene_size_dict[gene] = cnss.split(",")
    return gene_size_dict

def find(x,y):
    """returns ture is a set is not disjoint"""
    return x[1] >= y[0] or x[1] >= y[1]

def union(x,y):
    """joins two points if they are not disjoint"""
    if find(x,y):
        start = min(x[0],y[0])
        stop = max(x[1],y[1])
        return deque([(start,stop)])
    else:
        return deque([x,y])

def merged_feats(feats):
    prv_count = 1
    output = deque([])
    #print "feats2{0}".format(feats)
    while prv_count != len(feats) and len(feats) > 0:
        prv_count = len(feats)
        if len(feats) > 1:
            x = feats.popleft()
            y = feats.popleft()
            merged = union(x,y)
            output += merged
        else:
            if feats[0] not in output:
                output += feats
                feats.popleft()
       # print "output: {0}".format(output)
    return output

def recursive_merge(feats):
    new_feats = feats
    prv_count = 0
    while prv_count != len(new_feats) and len(new_feats) > 1:
        prv_count = len(new_feats)
        #print prv_count
        new_feats = merged_feats(new_feats)
    #print new_feats
    return new_feats

def recursive_merge_both(feats):
    evens = recursive_merge(feats)
    first_even = deque([evens[0]])
    evens.popleft()
    odds = recursive_merge(evens)
    if len(odds) > 1:
        first_odd = deque([odds[0]]) 
        odds.popleft()
        odds_two = recursive_merge(odds)
        list_both = list(first_even) + list(odds_two) + list(first_odd)
    else:
        list_both = list(first_even) + list(odds)
    list_both.sort()
    return list_both

def non_coding(feats,search_start,search_end):
    """find the inverse of genes ... the non_coding seq.. inverse the list"""
    non_cd_regions = []
    off_by_one = [(stop+1,(feats[i+1][0])) for
                    i,(start,stop) in enumerate(feats[:-1])]
    ###ii
    off_by_one.append((feats[-1][1],search_end))
    off_by_one.append((search_start,feats[0][0]))
    off_by_one.sort()
    return off_by_one

def random_noncoding(gene_names,bed,fasta_file,new_fasta_name):
    f = Fasta(fasta_file)
    new_fasta = open("{0}".format(new_fasta_name),"wb")
    for line in open(gene_names):
        gene = line.strip()
        gene = bed.accn(gene)
        search_start = gene["start"] - 12000
        search_end = gene["end"] + 12000

        feats = bed.get_features_in_region(gene["seqid"],search_start,search_end)
        #print search_start,search_end
        locs =[feat["locs"] for feat in feats]
        ##iilocs = [(feat['start'],feat['end']) for feat in feats]
        flatten_locs = chain.from_iterable(locs)
        flatten_locs = list(flatten_locs)
        ##iistart_stops =[(int(start),int(end)) for (start,end) in locs]
 
        start_stops =[(int(start),int(end)) for (start,end) in flatten_locs]
        start_stops_in_range = []
        for s_start,s_stop in start_stops:
            if s_start <= search_end and s_stop >= search_start:
                start_stops_in_range.append((s_start,s_stop))
        #print search_start,search_end
        #print start_stops_in_range
        start_stops_in_range.sort()
        #print start_stops
        if len(start_stops_in_range) > 1:
            start_stops_m = recursive_merge_both(deque(start_stops_in_range))
            #### need to merge when mask overlaps with gene
        else:
            start_stops_m = start_stops_in_range
        merged_feats = list(start_stops_m)
        merged_feats.sort()
        #print merged_feats
        off_by_one = non_coding(merged_feats,search_start,search_end)
        get_fasta(gene,off_by_one,f,new_fasta)
    new_fasta.close()

def get_fasta(accn,off_by_one,f,new_fasta):
    off_by_one.sort()
    seq = f[accn['seqid']]
    n = 0
    strand = accn['strand']
    for start,stop in off_by_one:
        n += 1 
        fasta = seq[start-1:stop-1]
        w = ">cns{1}_{0}\n".format(accn['accn'],n)
        if strand == '-':
            my_seq = fasta
            fasta = str(Seq(my_seq).reverse_complement())
        if len(fasta) == 0:
            #print start,stop,accn['accn']
            continue
        seq_w = "{0}\n".format(fasta)
        new_fasta.write(w)
        new_fasta.write(seq_w)

####### tair ##########
#x = random_noncoding('/Users/gt/Desktop/tmp.csv',Bed('/Users/gt/thaliana_v8.with_new_cns_mask.bed'),"/Users/gt/thaliana_v8.fasta","/Users/gt/thaliana_v8_control_SB.fasta")
x = random_noncoding('/Users/gt/Desktop/tmp.csv',Bed('/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_golden/thaliana_v8.with_new_cns_mask.bed'),"/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/thaliana_v8.fasta","/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_mine/thaliana_v8_control_SB.fasta")

######### rice,sorg,set #####
##### took out strand info used N to mask bed also ########


#x =
#random_noncoding('/Users/gt/Desktop/paper/G-box-seq/rice_rice/tmp.csv',Bed('/Users/gt/Desktop/paper/G-box-seq/rice.with_new_cns_mask.bed'),"/Users/gt/Desktop/paper/G-box-seq/rice_rice/rice_j.fasta","/Users/gt/Desktop/paper/G-box-seq/rice_rice/rice_rice_control_fasta")
#x = random_noncoding('/Users/gt/Desktop/tmp.csv',Bed('/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_mine/thaliana_v8.with_new_cns_mask.bed'),"/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/thaliana_v8.fasta","/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_mine/thaliana_v8_control_SB.fasta")
#x = random_noncoding('/Users/gt/Desktop/tmp.csv',Bed('/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_mine/thaliana_v8.with_new_cns_mask.bed'),"/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/thaliana_v8.fasta","/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_mine/thaliana_v8_control_SB.fasta")
#x = random_noncoding('/Users/gt/Desktop/tmp.csv',Bed('/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_mine/thaliana_v8.with_new_cns_mask.bed'),"/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/thaliana_v8.fasta","/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_mine/thaliana_v8_control_SB.fasta")
#x = random_noncoding('/Users/gt/Desktop/tmp.csv',Bed('/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_mine/thaliana_v8.with_new_cns_mask.bed'),"/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/thaliana_v8.fasta","/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_mine/thaliana_v8_control_SB.fasta")
#x = random_noncoding('/Users/gt/Desktop/tmp.csv',Bed('/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_mine/thaliana_v8.with_new_cns_mask.bed'),"/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/thaliana_v8.fasta","/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_mine/thaliana_v8_control_SB.fasta")
#






# x = random_noncoding('/Users/gt/Desktop/tmp.csv',Bed('/Users/gt/new/rice_j_setaria_n/rice_j_set.nolocaldups.with_new_cns_mask.bed'),"/Users/gt/data/paper4/rice_j.fasta","/Users/gt/Desktop/paper/rice_j_setaria_n_control_SB.fasta")
#x = random_noncoding('/Users/gt/Desktop/tmp.csv',Bed('/Users/gt/Desktop/footprint_cns/Tair_9_cns.with_new_cns.bed'),"/Users/gt/Desktop/footprint_cns/tair_9.fasta","/Users/gt/Desktop/footprint_cns/tair_9_control_nm.fasta")
###### rice_sorg #########
#dict_size = gene_size_dict('/Users/gt/Desktop/paper/tmp.tsv')
#x = random_noncoding(dict_size,Bed('/Users/gt/Desktop/paper/rice_j_sorg.nolocaldups.with_new_cns_mask.bed'))
#get_seq(x,"/Users/gt/data/paper4/rice_j.fasta","/Users/gt/Desktop/paper/rice_j_sorghum_n_control.fasta")


##### seq for cns
#handle = open("/Users/gturco/data/paper3/rice_b_sorghum_v1.cns.assigned_real.csv")
#fh = handle.read()
#cns_list = []
#for line in fh.split("\n")[:-1]:
#    if line[0] == "#": continue
#    cns_id,accn,seqid,start,end,strand = line.split(",")[:6]
#    cns_list.append((seqid,int(start),int(end)))
#
#len(cns_list)
#get_seq(cns_list,"/Users/gturco/data/paper3/rice_b.fasta","/Users/gturco/test_cns.fasta")
##
