import fisher 
import rpy2.robjects as ro
from collections import OrderedDict

def fisher_chi(study_c,study_p,control_c,control_p):
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

