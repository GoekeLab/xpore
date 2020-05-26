import numpy as np
from numpy import newaxis
import scipy.special
import scipy.stats


class GMM(object):  
    """
    1D multi-sample 2-Gaussian mixture model.
    """
    def __init__(self, method=None, data={'x': None, 'y': None, 'condition_names': None, 'run_names': None}, inits={'info': None, 'nodes': {'x': None, 'y': None, 'w': None, 'mu_tau': None, 'z': None}}, priors={'mu_tau': None, 'w': None},kmer_signal=None):
        """
        Parameters
        ----------
        method: dict
            Config file
        data: dict
            Sth.
        inits: dict
            Sth.
        priors:
            Sth.
        """
        
        self.nodes = dict()
        self.aux_params = dict()
        self.K = 2  # modified and unmodified
        self.method = method
        self.kmer_signal = kmer_signal
        
        self.__init_info(inits['info'])

        data_node_x = None
        if (data['x'] is not None) and (data['y'] is not None) and (data['r'] is not None):
            
            self.info['n_reads'] = len(data['y'])

            # inits
            location_0 = priors['mu_tau']['location'][0]
            if location_0 < np.percentile(data['y'], q=50):
                location_1 = np.percentile(data['y'], q=75)
            else:
                location_1 = np.percentile(data['y'], q=25)
            locations = np.array([location_0, location_1])
            # print(locations)

            # kmeans = KMeans(n_clusters=self.K).fit(data['y'][:,newaxis])
            # locations = kmeans.cluster_centers_.flatten()

            inits['nodes']['mu_tau'] = {'location': locations, 'lambda': priors['mu_tau']['lambda'], 'alpha': priors['mu_tau']['alpha'], 'beta': priors['mu_tau']['beta']}
            
            if self.method['pooling']:
                inits['nodes']['x'] = {'group_names': data['condition_names']}
                data_node_x = data['x']
            else:
                inits['nodes']['x'] = {'group_names': data['run_names']}
                data_node_x = data['r']

            # additional nodes for postprocessing analysis
            self.nodes['y_condition_names'] = Constant(data=data['y_condition_names'])
            self.nodes['y_run_names'] = Constant(data=data['y_run_names'])

        # Define the graph
        self.nodes['x'] = Constant(data=data_node_x, inits=inits['nodes']['x'])   # NG
        self.n_groups = len(self.nodes['x'].params['group_names'])
        self.nodes['w'] = Dirichlet(dim=(self.n_groups, self.K), inits=inits['nodes']['w'], priors=priors['w'])
        self.nodes['mu_tau'] = UnivariateNormalGamma(dim=(self.K), inits=inits['nodes']['mu_tau'], priors=priors['mu_tau'])
        self.nodes['z'] = Bernoulli(dim=(self.info['n_reads'], self.K), parents={'w': self.nodes['w'], 'x': self.nodes['x']}, inits=inits['nodes']['z'])
        self.nodes['y'] = UnivariateNormalMixture(parents={'z': self.nodes['z'], 'mu_tau': self.nodes['mu_tau'], 'x': self.nodes['x']}, data=data['y'], inits=inits['nodes']['y'])

    def __init_info(self, inits):
        if inits is None:
            self.info = {'n_reads': 0, 'log_elbos': [], 'converged': False, 'n_iterations': -1, 'convergence_ratio': -1}
        else:
            self.info = inits

    def __compute_log_elbo(self):
        """
        Compute log ELBO.
        """
        log_elbo = self.nodes['y']._log_likelihood()
        log_elbo += self.nodes['z']._log_prob_prior()
        log_elbo -= self.nodes['z']._log_prob_posterior()
        log_elbo += self.nodes['w']._log_prob_prior()
        log_elbo -= self.nodes['w']._log_prob_posterior()
        log_elbo += self.nodes['mu_tau']._log_prob_prior(var='normal')
        log_elbo += self.nodes['mu_tau']._log_prob_prior(var='gamma')
        log_elbo -= self.nodes['mu_tau']._log_prob_posterior(var='normal')
        log_elbo -= self.nodes['mu_tau']._log_prob_posterior(var='gamma')

        return log_elbo

    def fit(self):
        """
        Fit.
        """
        error = False
        converged = False
        # self.__init_params()
        log_elbo_old = self.__compute_log_elbo()  # -np.inf
        self.nodes['y']._update()

        for iteration in range(1, self.method['max_iters']):
            self.nodes['z']._update(children={'y': self.nodes['y']})  # every time z is updated, y needs to be done too. Update z so that z is not randomly drawn.
            self.nodes['y']._update()
            self.nodes['w']._update(children=self.nodes['y'].params)
            self.nodes['mu_tau']._update(children={**self.nodes['y'].params, **{'x': self.nodes['x']}})

            log_elbo_new = self.__compute_log_elbo()
            error = np.round(log_elbo_new - log_elbo_old, 8) < 0
            if error:
                print('ERROR: log ELBO is decreasing ...')
                print(iteration, log_elbo_old, log_elbo_new, error)
                error = True
                break
            else:
                self.info['log_elbos'] += [log_elbo_old]
                diff = (log_elbo_new - log_elbo_old)/np.abs(log_elbo_old)
                log_elbo_old = log_elbo_new
                # print('Iteration %d:' %(iteration),'Convergence ratio =',diff)
                if (diff < self.method['stopping_criteria']):
                    converged = True
                    break
        self.info['n_iterations'] = iteration
        self.info['converged'] = converged
        self.info['convergence_ratio'] = diff
        return self

################################
##### Model for each node ######
################################


class Constant(object):
    """
    Constant node.
    """
    def __init__(self, data=None, inits=None):
        self.params = inits
        self.data = data


class UnivariateNormalMixture(object):
    def __init__(self, parents=None, data=None, inits=None):

        self.parents = parents
        self.data = data

        self.params = dict.fromkeys(['N', 'mean', 'variance'])
        self.__initialise_params(inits)

    def __initialise_params(self, inits):
        if inits is not None:
            self.params = inits

    def _log_likelihood(self):
        n_reads = len(self.data)
        mu_tau_nomal_means = np.repeat(self.parents['mu_tau'].expected(var='normal')[newaxis, :], repeats=n_reads, axis=-2)  # N,K
        mu_tau_normal_variances = np.repeat(self.parents['mu_tau'].variance(var='normal')[newaxis, :], repeats=n_reads, axis=-2)  # N,K
        mu_tau_gamma_means = np.repeat(self.parents['mu_tau'].expected(var='gamma')[newaxis, :], repeats=n_reads, axis=-2)  # N,K
        mu_tau_log_gamma = np.repeat(self.parents['mu_tau'].expected(var='gamma', inner_func='log')[newaxis, :], repeats=n_reads, axis=-2)  # N,K
        y = self.data[..., newaxis]  # N,1

        log_likelihood = self.parents['z'].expected() * 0.5 * (mu_tau_log_gamma - np.log(2*np.pi) - mu_tau_gamma_means * (y**2 - 2*y * mu_tau_nomal_means + mu_tau_nomal_means**2 + mu_tau_normal_variances))  # N,K
        log_likelihood *= (np.sum(self.parents['x'].data, axis=-1)[..., newaxis] > 0)  # N,K
        return np.sum(np.sum(log_likelihood, axis=-1), axis=-1)  # sum across components and reads respectively => 1

    def expected(self, inner_func=None):
        return self.params['mean']

    def variance(self, inner_func=None):
        return self.params['variance']

    def N(self, inner_func=None):
        return self.params['N']

    def _update(self):
        n_reads = len(self.data)
        self.params['N'] = np.sum(self.parents['x'].data[..., newaxis] * self.parents['z'].expected()[:, newaxis, :], axis=-3)  # sum across reads => G,K
        N = np.sum(self.params['N'], axis=-2)  # sum across groups => K
        N[N < 1e-30] = 0.  # Hack: Avoid overflow in later calculation steps.

        N_inverse = np.divide(1., N, out=np.zeros_like(N), where=N != 0.)
        y = self.data[..., newaxis]  # N,1
        self.params['mean'] = N_inverse * np.sum(self.parents['z'].expected()*y, axis=-2)  # sum across reads => K

        residuals = (np.sum(self.parents['x'].data, axis=-1) > 0)[..., newaxis]*(y - np.repeat(self.params['mean'][newaxis, :], repeats=n_reads, axis=-2))  # N,K

        self.params['variance'] = N_inverse * np.sum(self.parents['z'].expected()*(residuals**2), axis=-2)  # sum across reads => K


class Bernoulli(object):
    def __init__(self, dim, parents=None, data=None, inits=None):
        self.params = dict()
        self.params['prob'] = np.full(dim, np.nan)
        self.params['ln_prob'] = np.full(dim, np.nan)
        self.__initialise_params(inits)

        self.parents = parents
        self.data = data

    def __initialise_params(self, inits=None):
        if inits is None:
            # self.params['prob'][:] = 0.5
            # self.params['ln_prob'][:] = np.log(0.5)

            self.params['prob'] = np.random.random(self.params['prob'].shape)
            self.params['prob'] /= np.sum(self.params['prob'], axis=1)[:, newaxis]
            self.params['ln_prob'] = np.log(self.params['prob'])
        else:
            self.params = inits

    def _log_prob_prior(self):
        res = self.expected() * np.sum(self.parents['x'].data[..., newaxis] * self.parents['w'].expected(inner_func='log')[newaxis, :, :], axis=-2)  # sum across groups => N,K
        return np.sum(np.sum(res, axis=-1), axis=-1)  # 1

    def _log_prob_posterior(self):
        res = self.expected() * np.sum(self.parents['x'].data[:, newaxis, :] * self.params['ln_prob'][..., newaxis], axis=-1)  # N,K
        return np.sum(np.sum(res, axis=-1), axis=-1)  # 1

    def expected(self, inner_func=None):
        return self.params['prob']

    # def variance(self):
    #     return self.params['prob']*(1-self.params['prob']) #scipy.stats.bernoulli.var(self.params['concentration']) #

    def _update(self, children):
        # print('#'*5,'Updating z','#'*5)
        n_reads = len(children['y'].data)
        mu_tau_nomal_means = np.repeat(children['y'].parents['mu_tau'].expected(var='normal')[newaxis, :], repeats=n_reads, axis=-2)  # N,K
        mu_tau_normal_variances = np.repeat(children['y'].parents['mu_tau'].variance(var='normal')[newaxis, :], repeats=n_reads, axis=-2)  # N,K
        mu_tau_gamma_means = np.repeat(children['y'].parents['mu_tau'].expected(var='gamma')[newaxis, :], repeats=n_reads, axis=-2)  # N,K
        mu_tau_expected_log_gamma = np.repeat(children['y'].parents['mu_tau'].expected(var='gamma', inner_func='log')[newaxis, :], repeats=n_reads, axis=-2)  # N,K

        y = children['y'].data[..., newaxis]  # N,1
        ln_rho = 0.5*(mu_tau_expected_log_gamma - np.log(2*np.pi) - mu_tau_gamma_means * (y**2 - 2*y * mu_tau_nomal_means + mu_tau_nomal_means**2 + mu_tau_normal_variances))

        ln_rho += np.sum(self.parents['x'].data[..., newaxis]*self.parents['w'].expected(inner_func='log')[newaxis, :, :], axis=-2)  # sum across groups => N,K

        ln_rho_diff = ln_rho[:, 1]-ln_rho[:, 0]
        prob = 1./(1.+np.exp(np.clip(ln_rho_diff, -500, 500)))  # N,1 => N
        self.params['prob'] = np.stack([prob, 1-prob], axis=-1)
        self.params['ln_prob'] = ln_rho - scipy.special.logsumexp(ln_rho, axis=-1)[..., newaxis]  # N,K


class Dirichlet(object):
    def __init__(self, dim, parents=None, data=None, inits=None, priors=None):  # dim - [,n_categories]

        self.priors = dict()
        self.priors['concentration'] = np.full(dim, np.nan)
        self.__set_priors(priors)

        self.params = dict()
        self.params['concentration'] = np.full(dim, np.nan)
        self.__initialise_params(inits)

        self.parents = parents
        self.data = data

    def __initialise_params(self, inits=None):
        if inits is None:
            self.params['concentration'][:] = self.priors['concentration'][:]
        else:
            self.params = inits

    def __set_priors(self, priors=None):
        if priors is None:
            self.priors['concentration'][:] = 1
        else:
            self.priors = priors

    def _log_prob_prior(self):
        res = Dirichlet.__log_C(self.priors['concentration']) + np.sum((self.priors['concentration']-1)*self.expected(inner_func='log'), axis=-1)  # sum k => G
        return np.sum(res, axis=-1)  # 1

    def _log_prob_posterior(self):
        res = Dirichlet.__log_C(self.params['concentration']) + np.sum((self.params['concentration']-1)*self.expected(inner_func='log'), axis=-1)  # sum k => G
        return np.sum(res, axis=-1)  # 1

    def expected(self, inner_func=None):
        if inner_func == 'log':
            res = scipy.special.digamma(self.params['concentration'])  # G,K
            res -= scipy.special.digamma(np.sum(self.params['concentration'], axis=-1))[..., newaxis]  # G,K / G,1 => G,K
            return res
        else:
            return self.params['concentration']/np.sum(self.params['concentration'], axis=-1)[..., newaxis]  # G,K

    def variance(self):  # To check, still wrong but unused.
        return (self.expected()*(1-self.expected())) / (np.sum(self.params['concentration'], axis=-1)[..., newaxis]+1)

    def _update(self, children):
        # print('#'*5,'Updating w','#'*5)
        self.params['concentration'] = self.priors['concentration'] + children['N']

    def __log_C(alpha):
        return scipy.special.gammaln(np.sum(alpha, axis=-1)) - np.sum(scipy.special.gammaln(alpha), axis=-1)


class UnivariateNormalGamma(object):
    def __init__(self, dim, parents=None, data=None, inits=None, priors=None):

        self.priors = dict.fromkeys(['location', 'lambda', 'alpha', 'beta'], np.full(dim, np.nan))  # alpha = shape, beta=rate
        self.__set_priors(priors)

        self.params = dict.fromkeys(['location', 'lambda', 'alpha', 'beta'], np.full(dim, np.nan))  # alpha = shape, beta=rate
        self.__initialise_params(inits)
        self.parents = parents
        self.data = data

    def __initialise_params(self, inits):
        if inits is None:
            for param in self.params:
                self.params[param] = self.priors[param][:]
            self.params['location'][:] = np.random.normal(loc=self.priors['location'][:], scale=1)

        else:
            self.params = inits

    def __set_priors(self, priors=None):
        # if priors is None:
        #     self.priors['location'][:] = 0. # Todo: set to mean
        #     self.priors['lambda'][:] = 1.
        #     self.priors['alpha'][:] = 1.
        #     self.priors['beta'][:] = 1.
        # else:
        self.priors = priors

    def _log_prob_prior(self, var='normal'):
        if var == 'normal':
            res = 0.5 * (self.expected(var='gamma', inner_func='log') + np.log(self.priors['lambda']) - np.log(2*np.pi) - self.priors['lambda']*self.expected(var='gamma') * (self.expected(var='normal')**2 + self.variance(var='normal') - 2*self.expected(var='normal')*self.priors['location'] + self.priors['location']**2))  # K

        else:  # gamma
            res = self.priors['alpha']*np.log(self.priors['beta'])-scipy.special.gammaln(self.priors['alpha'])+(self.priors['alpha']-1)*self.expected(var='gamma', inner_func='log') - self.priors['beta']*self.expected(var='gamma')  # K

        return np.sum(res, axis=-1)  # 1

    def _log_prob_posterior(self, var='normal'):
        if var == 'normal':
            res = 0.5 * (self.expected(var='gamma', inner_func='log') + np.log(self.params['lambda']) - np.log(2*np.pi) - self.params['lambda']*self.expected(var='gamma') * (self.expected(var='normal')**2 + self.variance(var='normal') - 2*self.expected(var='normal')*self.params['location'] + self.params['location']**2))  # K

        else:  # gamma
            res = self.params['alpha']*np.log(self.params['beta'])-scipy.special.gammaln(self.params['alpha'])+(self.params['alpha']-1)*self.expected(var='gamma', inner_func='log') - self.params['beta']*self.expected(var='gamma')  # K

        return np.sum(res, axis=-1)  # 1

    def expected(self, var='normal', inner_func=None):
        if inner_func == 'log':
            if var == 'gamma':  # tau
                return scipy.special.digamma(self.params['alpha']) - np.log(self.params['beta'])  # K

        else:
            if var == 'normal':
                return self.params['location']
            else:  # gamma
                return self.params['alpha'] / self.params['beta']

    def variance(self, var='normal'):
        if var == 'normal':
            return 1./(self.params['lambda']*self.expected(var='gamma'))
        # else: # gamma
        #     return self.params['alpha']/(self.params['beta']**2)

    def _update(self, children):
        # print('#'*5,'Updating mu,tau','#'*5)
        N = np.sum(children['N'], axis=-2)  # sum across groups => K
        self.params['lambda'] = self.priors['lambda'] + N  # K
        self.params['location'] = (1./self.params['lambda']) * (N*children['mean'] + self.priors['lambda']*self.priors['location'])  # K
        self.params['alpha'] = self.priors['alpha'] + 0.5*N  # K
        self.params['beta'] = self.priors['beta'] + 0.5*(N*children['variance'] + ((self.priors['lambda']*N) / (self.priors['lambda']+N))*(children['mean']-self.priors['location'])**2)  # K
