library('trendsceek')

pp = pos2pp(pos)
log.fcn = log10
pp = set_marks(pp, t(10^value), log.fcn = log.fcn)

pp2plot = pp_select(pp)

##set parameters
nrand = 100
ncores = 1
##run
trendstat_list = trendsceek_test(pp2plot, nrand, ncores)

head(trendstat_list[['supstats_wide']])
alpha = 0.1 ##Benjamini-Hochberg
sig_list = extract_sig_genes(trendstat_list, alpha)
sig_list
