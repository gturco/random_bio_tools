import sys
from pyfasta import Fasta
from Bio.SeqUtils import GC
import re
#x = 'AAAGC'
#40.0=GC(x)
def genome_contenct_stats(fasta_path):
    f = Fasta(fasta_path)
    g_box_total = []
    for seqid in f.keys():
        seq = f[seqid][:]
        g_boxs = len(re.findall('CACGTG',seq,flags=re.IGNORECASE))
        g_box_total.append(g_boxs)
    print >>sys.stderr, "total gboxes:{0}".format(sum(g_box_total))

def gc_stats(fasta_file,header,seven_mer_seq):
    """takes a 'header' seq fasta file and uses headers as names
    returns a dictionary with gc cont, g_box number, masked regions for each
    header/seq"""
    fh = open(fasta_file)
    fasta = fh.read()
    gene_content = {}
    fasta_lines = fasta.split('\n')[:-1]
    seq_len = []
    for pos,line in enumerate(fasta_lines):
        if re.match('^{0}'.format(header),line):
            name = line[1:]
            seq = fasta_lines[pos+1]
            gene_content[name]={}
            gc = GC(seq)
            gene_content[name]['gc']=gc
            seven_mer = len(re.findall(seven_mer_seq,seq,flags=re.IGNORECASE))
            g_boxs = len(re.findall('CACGTG',seq,flags=re.IGNORECASE))
            gene_content[name]['g_boxs'] = g_boxs
            gene_content[name]['seven_mer'] = seven_mer
            masked = len(re.findall('X|N',seq,flags=re.IGNORECASE))
            if len(seq) == 0:
                #print pos
                gene_content[name]['masked'] = 100
            else:
                gene_content[name]['masked'] = float(masked/len(seq)*100)
            seq_len.append(len(seq)-masked)
    print "seq len: {0}".format(sum(seq_len))
    return gene_content, sum(seq_len)

def stats(fasta_file,header,seven_mer_seq):
    stat_dict, seq_len = gc_stats(fasta_file, header,seven_mer_seq)
    # for percent GC need to remove all genes that are masked
    gc_percent = []
    g_box_count = []
    seven_mer_count = []
    for key in stat_dict.keys():
        if stat_dict[key]['masked'] > 0: continue
        g_box_count.append(stat_dict[key]['g_boxs'])
        seven_mer_count.append(stat_dict[key]['seven_mer'])
        # for percent GC need to remove all genes that are masked
        gc_percent.append(stat_dict[key]['gc'])
    seven_mer_freq = sum(seven_mer_count)
    g_box_freq = sum(g_box_count)
    gc_percent_all = float(sum(gc_percent))/len(gc_percent)
    print >>sys.stderr, "g_boxs:{0} \t gc_percent:{1} \t number_of_genes {2} \t seven_mer {3}".format(g_box_freq,gc_percent_all, len(gc_percent),seven_mer_freq)
    return seven_mer_freq, seq_len

#genome_contenct_stats("/Users/gturco/data/paper3/rice_b.fasta")
stats("/Users/gturco/test.fasta",">c","ACTGGGGG")

#g_box_all = []
#seq_len = []
#for fasta in range(1,7):
#    handle = open("/Users/gturco/{0}.fasta".format(fasta))
#    fh = handle.read()
#    seq = fh.split("\n")[0]
#    g_boxs = len(re.findall('CACGTG',seq,flags=re.IGNORECASE))
#    g_box_all.append(g_boxs)
#    #print g_boxs, fasta
#    masked = len(re.findall('X|N',seq,flags=re.IGNORECASE))
#    seq_len.append(len(seq)-masked)
#print "{0}\n{1}".format(sum(g_box_all), sum(seq_len))
