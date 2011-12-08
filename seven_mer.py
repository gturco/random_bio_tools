from fisher_chisq import fisher_chi
from paper import stats
from itertools import product

x = list(product(['A','C','T','G'], repeat = 7))


def main(study_fasta_file,control_fast_file,out_file, len_mer):
    handle = open(out_file,"wb")
    record = []
    seven_mers = list(product(['A','C','T','G'], repeat = len_mer))

    for seven_mer in seven_mers:
        seq = "".join(seven_mer)
        study_c,study_p = stats(study_fasta_file,">cns",seq)
        control_c,control_p = stats(control_fast_file,">cns",seq)
        res, enrich = fisher_chi(study_c,study_p,control_c,control_p)
        print res
        record.append((res,enrich,seven_mer,study_c,study_p,control_c,control_p))
    record.sort()
    for res,enrich,seven_mer,study_c,study_p,control_c,control_p in record:
        w = "{0}\t{1}\t{2}/{3}\t{4}/{5}\t{6}\n".format(seven_mer,enrich,study_c,study_p,control_c,control_p,res)
        handle.write(w)


main("/Users/gturco/data/paper3/rice_b_setaria64.cns_real.fasta","/Users/gturco/data/paper3/rice_b_setaria64.cns_random.fasta","/Users/gturco/test_mers",7)
