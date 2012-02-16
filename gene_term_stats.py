import fisher 
import rpy2.robjects as ro
from collections import Counter

def fisher_chi(study_c,study_p,control_c,control_p):
    """takes study count, population count, control count and control total"""
    if study_p <= 1000 and control_p <= 1000:
        p = fisher.pvalue(study_c,study_p,control_c,control_p)
        res = [p.two_tail]
    else:
        v = ro.IntVector([study_c, study_p, control_c, control_p])
        m = ro.r['matrix'](v, 2, 2)
        res = ro.r['chisq.test'](m)[2]
    enrichment = 'e' if 1.0* study_c/study_p > 1.0 * control_c / control_p else 'p'
    return res[0],enrichment

#res= fisher_chi(100,200,30,400)
##
def count_terms(term_file):
    """ input file were each line is a term for either pop or study.
    output term count along with total number of terms in file"""
    cnt = Counter()
    total = 0
    for line in open():
        word = line.strip()
        cnt[word] += 1
        total += 1
    return cnt,total

def term_stats(study, population,out_fh):
    study_count,study_all = count_terms(study)
    pop_count,pop_all = count_terms(population)
    write_file = open(out_fh,"wb")
    for term in study_count.keys():
        study_c = study_count[term]
        try:
            pop_c = pop_count[term]
        except KeyError:
            raise "population file does not have all terms of study!!!"
        study_n = study_all - study_c
        pop_n = (pop_all - pop_c - study_all) + study_c
        res,enrichment = fisher_chi(study_c,study_n,pop_c,pop_n)
        write_file.write("{0}\t{1}/{2}\t{3}/{4}\t{5}\t{6}\n".format(term,study_c,study_all,pop_c,pop_all,res,enrichment))
    write_file.close()

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("-p", dest="pop", help="file contains list of terms for each gene in  popultation")
    parser.add_option("-s", dest="study", help="file contains list of terms for each gene in study")
    parser.add_option("-o", dest="outfile",help="file that stats are wrote to")
    (options, _) = parser.parse_args()

    term_stats(options.pop,options.study,options.outfile)
