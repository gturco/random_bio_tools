import sys
from pyfasta import Fasta
from Bio.SeqUtils import Seq
import re
from collections import defaultdict
from itertools import product
#x = 'AAAGC'
#40.0=GC(x)


def get_gbox(my_seq,m):
    g_boxs = len(re.findall(m,my_seq,flags=re.IGNORECASE))
    return g_boxs

def get_neg(my_seq,m):
    fasta = str(Seq(my_seq).reverse_complement())
    g_boxs = get_gbox(fasta,m)
    return g_boxs

def motif_count(fasta_file,motif,expand,out):
    seq_len = []
    pos = defaultdict(list)
    neg = defaultdict(list)
    fastah = open(fasta_file)
    fasta = fastah.read()
    fasta_lines = fasta.split('\n')
    if fasta_lines[-1] == '':
        fasta_lines = fasta_lines[:-1]
    for i,line in enumerate(fasta_lines[:-1]):
        if '>'in line:
            name = '{0}_{1}'.format(line[1:],i)
            seq = fasta_lines[i+1]
            seq_len.append(len(seq))
            motifs = generate_ends(motif,expand) + [motif]
            for m in motifs:
                pos[m].append(get_gbox(seq,m))
                neg[m].append(get_neg(seq,m))
    return pos,neg,sum(seq_len)

def main(control,test,motif,expand,outfile):
    out = open(outfile,'wb')
    c_pos,c_neg,c_seq = motif_count(control,motif,expand,out)
    t_pos,t_neg,t_seq = motif_count(test,motif,expand,out)
    write_results(c_pos,c_neg,c_seq,t_pos,t_neg,t_seq,out)


def write_results(c_pos,c_neg,c_seq,t_pos,t_neg,t_seq,out):
    for key in t_pos.keys():
        a = ',{0}/{1},{2}'.format(sum(t_pos[key]),t_seq,sum(t_pos[key])/float(t_seq))
        b = '{0}/{1},{2}'.format(sum(t_neg[key]),t_seq,sum(t_neg[key])/float(t_seq))
        c = '{0}/{1},{2}'.format(sum(c_pos[key]),c_seq,sum(c_pos[key])/float(c_seq))
        d = '{0}/{1},{2},'.format(sum(c_neg[key]),c_seq,sum(c_neg[key])/float(c_seq))
        tpos_avg = (sum(t_pos[key])/float(t_seq))
        cpos_avg = (sum(c_pos[key])/float(c_seq))
        tneg_avg = (sum(t_neg[key])/float(t_seq))
        cneg_avg = (sum(c_neg[key])/float(c_seq))
        e = '{0},{1}\n'.format(tpos_avg/cpos_avg,tneg_avg/cneg_avg)
        print tpos_avg/cpos_avg , key
        line = key + a + b + c + d + e
       
        out.write(line)
    out.close()



def generate_ends(motif,ends):
    palindromes = []
    for end in product('ACGT',repeat=ends):
        seq = ''.join(end)
        rc = str(Seq(seq).reverse_complement())
        palindromes.append(seq + motif + rc)
    return palindromes



#main('control_fasta/rice_rice/all.fasta','/Users/gt/Desktop/paper/G-box-seq/rice_rice/rice_j_rice_j.cns_real.fasta','CACGTG',2,'out')
#main('control_fasta/tair_8/all.fasta','/Users/gt/Desktop/freelinglab/genomes/tair 8/tair_8/tair_8_golden/testing','CACGTG',2,'out')
#main('control_fasta/sorg_sorg/all.fasta','/Users/gt/Desktop/paper/G-box-seq/sorg_sorg/sorghum_n_sorghum_n.cns_real.fasta','CACGTG',2,'out_sorg')
#main('control_fasta/set_set/all.fasta','/Users/gt/Desktop/paper/G-box-seq/set_set/setaria_n_setaria_n.cns_real.fasta','CACGTG',2,'out_set')
#main('control_fasta/rice_set/all.fasta','/Users/gt/Desktop/paper/G-box-seq/rice_set/rice_only.fasta','CACGTG',2,'out_rice_set')
main('control_fasta/rice_sorg/all.fasta','/Users/gt/Desktop/paper/G-box-seq/rice_sorg/rice_j_only_cns.fasta','CACGTG',2,'out_rice_sorg')




#
