import scipy.stats
import numpy as np

class StatsTest(object):
    def __init__(self,data):
        if self.__isok(data):
            self.data = [data['y'][data['x'][:,0]==1],data['y'][data['x'][:,1]==1]] 
        else:
            self.data = None

    def fit(self,method):
        pval = np.nan
        if self.data is not None:
            if method == 't-test':
                stats,_ = scipy.stats.ttest_ind(self.data[0],self.data[1])
                df = len(self.data[0]) + len(self.data[1]) - 2
                pval = scipy.stats.t.sf(np.abs(stats), df)*2 # also defined as 1 - cdf, but sf is sometimes more accurate. #(1 - scipy.stats.t.cdf(abs(stat), df)) * 
        return pval
    
    def __isok(self,data):
        return data['x'].shape[1] == 2 # Check if only two conditions.


