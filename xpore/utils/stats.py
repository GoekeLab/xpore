import scipy.stats
from math import fabs, sqrt, log, erf
import numpy as np
# from ..exts.statistics import NormalDist


def z_test(y1, y2, n1, n2):
    # two-tailed
    p1 = y1.mean()
    p2 = y2.mean()
    n1 = n1.mean().round()
    n2 = n2.mean().round()
    se = np.sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2)
    z = (p1 - p2) / se
    return z, scipy.stats.norm.sf(abs(z))*2


def calc_prob_overlapping(means, variances):
    sigma1, sigma2 = np.sqrt(variances)
    mu1, mu2 = means
    return NormalDist(mu=mu1, sigma=sigma1).overlap(NormalDist(mu=mu2, sigma=sigma2))

class NormalDist(object):
    "Modified from the library statistics.py (Authors)(https://github.com/python/cpython/blob/3.8/Lib/statistics.py) which is under GPL-compatible license."
    def __init__(self, mu=0.0, sigma=1.0):
        self._mu = mu
        self._sigma = sigma
    
    def cdf(self, x):
        return 0.5 * (1.0 + erf((x - self._mu) / (self._sigma * sqrt(2.0))))

    def overlap(self, another):
        """
        Compute the overlapping coefficient (OVL) between two normal distributions.
        Measures the agreement between two normal probability distributions.
        Returns a value between 0.0 and 1.0 giving the overlapping area in
        the two underlying probability density functions.
        """ 
        # See: "The overlapping coefficient as a measure of agreement between
        # probability distributions and point estimation of the overlap of two
        # normal densities" -- Henry F. Inman and Edwin L. Bradley Jr
        # http://dx.doi.org/10.1080/03610928908830127
        X, Y = self, another
        if (Y._sigma, Y._mu) < (X._sigma, X._mu):  # sort to assure commutativity
            X, Y = Y, X
        X_var, Y_var = X.variance, Y.variance
        
        dv = Y_var - X_var
        dm = fabs(Y._mu - X._mu)
        if not dv:
            return 1.0 - erf(dm / (2.0 * X._sigma * sqrt(2.0)))
        a = X._mu * Y_var - Y._mu * X_var
        b = X._sigma * Y._sigma * sqrt(dm**2.0 + dv * log(Y_var / X_var))
        x1 = (a + b) / dv # intersection point 1
        x2 = (a - b) / dv # intersection point 2
        
        p_overlap = 1.0 - (fabs(Y.cdf(x1) - X.cdf(x1)) + fabs(Y.cdf(x2) - X.cdf(x2)))
        list_cdf_at_intersections = [X.cdf(x1),Y.cdf(x1),X.cdf(x2),Y.cdf(x2)]
        
        return p_overlap,list_cdf_at_intersections

    @property
    def mean(self):
        "Arithmetic mean of the normal distribution."
        return self._mu

    @property
    def stdev(self):
        "Standard deviation of the normal distribution."
        return self._sigma

    @property
    def variance(self):
        "Square of the standard deviation."
        return self._sigma ** 2.0

