.. _outputtable:

Output table description
=========================

==========================================  ========================================================================================================================================
Column name                                 Description
==========================================  ========================================================================================================================================
id                                          transcript or gene id
position                                    transcript or gene position
kmer                                        5-mer where modified base sits in the middle if modified
diff_mod_rate_<condition1>_vs_<condition2>  differential modification rate between condition1 and condition2 (modification rate of condition1 - modification rate of condition2)
z_score_<condition1>_vs_<condition2>        z score obtained from z-test of the differential modification rate
pval_<condition1>_vs_<condition2>           significance level from z-test of the differential modification rate
mod_rate_<condition>-<replicate>            modification rate of a replicate in the condition
mu_unmod                                    inferred mean of the unmodified RNAs distribution
mu_mod                                      inferred mean of the modified RNAs distribution
sigma2_unmod                                inferred sigma^2 of the unmodified RNAs distribution
sigma2_mod                                  inferred sigma^2 of the modified RNAs distribution
conf_mu_unmod                               confidence level of mu_unmod compared to the unmodified reference signal
conf_mu_mod                                 confidence level of mu_unmod compared to the unmodified reference signal
mod_assignment                              lower if mu_mod < mu_unmod and higher if mu_mod > mu_unmod
==========================================  ========================================================================================================================================
