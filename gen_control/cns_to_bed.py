from flatfeature import Bed
from make_genelist import merge_flat
from pyfasta import Fasta
import re

def cns_to_bed(cns_list,cns_bed_file):
    cns_bed = open("{0}".format(cns_bed_file),"wb")
    handle = open(cns_list)
    fh = handle.read()
    for line in fh.split("\n")[:-1]:
        if line[0]== "#": continue
        #print line.split(",")[:5]
        cns_id,accn,seqid,start,end,strand = line.split(",")[:6]
        w = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t.\t.\t.\t1\t{6}\t0\n'.format(seqid,start,end,cns_id,(int(end)-int(start)),strand,(int(end)-int(start)+1))
        cns_bed.write(w)


#cns_to_bed("/Users/gt/data/paper4/rice_j_setaria_n/rice_j_setaria_n.cns.assigned_real.csv","/Users/gt/data/paper4/rice_j_setaria_n/rice_j_setaria_n_cns.bed")
#merge_flat("/Users/gt/data/paper4/rice_j_setaria_n/rice_j_set.nolocaldups.with_new_cns.bed",Bed("/Users/gt/data/paper4/rice_j_setaria_n/rice_j.nolocaldups.with_new.all.local"),Bed("/Users/gt/data/paper4/rice_j_setaria_n/rice_j_setaria_n_cns.bed"))
cns_to_bed("/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorghum_n.cns.assigned_real.csv","/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorghum_n_cns.bed")
merge_flat("/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorg.nolocaldups.with_new_cns.bed",Bed("/Users/gt/data/paper4/rice_j_sorghum_n/rice_j.nolocaldups.with_new.all.local"),Bed("/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorghum_n_cns.bed"))


def mask_to_bed(fasta_file, mask_bed_name):
    "creates a bed file of the start and stops of masked seqs"
    mask_bed = open(mask_bed_name,"wb")
    f= Fasta(fasta_file)
    mask_id = 1
    for seqid in f.keys():
        seq = f[seqid][:]
        for m in re.finditer("X+",seq):
            mask_id = mask_id + 1
            w = '{0}\t{1}\t{2}\t{3}\t{4}\t+\t.\t.\t.\t1\t{5}\t0\n'.format(seqid,m.start(),m.end(),"mask_id {0}".format(mask_id),(m.end()-m.start()),(m.end()-m.start()+1))
            mask_bed.write(w)

#mask_to_bed("/Users/gt/data/paper4/rice_j.fasta","/Users/gt/data/paper4/rice_j_setaria_n/rice_j_set.masked_bed")
#merge_flat("/Users/gt/data/paper4/rice_j_setaria_n/rice_j_set.nolocaldups.with_new_cns_mask.bed",Bed("/Users/gt/data/paper4/rice_j_setaria_n/rice_j_set.nolocaldups.with_new_cns.bed"),Bed("/Users/gt/data/paper4/rice_j_setaria_n/rice_j_set.masked_bed"))
mask_to_bed("/Users/gt/data/paper4/rice_j.fasta","/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorg.masked_bed")
merge_flat("/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorg.nolocaldups.with_new_cns_mask.bed",Bed("/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorg.nolocaldups.with_new_cns.bed"),Bed("/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorg.masked_bed"))

