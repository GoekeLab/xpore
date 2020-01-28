import scipy.stats
import numpy
from ..exts.statistics import NormalDist


def z_test(y1, y2, n1, n2):
    # two-tailed
    p1 = y1.mean()
    p2 = y2.mean()
    n1 = n1.mean().round()
    n2 = n2.mean().round()
    se = numpy.sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2)
    z = (p1 - p2) / se
    return z, scipy.stats.norm.sf(abs(z))*2


def calc_prob_overlapping(means, variances):
    sigma1, sigma2 = numpy.sqrt(variances)
    mu1, mu2 = means
    return NormalDist(mu=mu1, sigma=sigma1).overlap(NormalDist(mu=mu2, sigma=sigma2))
