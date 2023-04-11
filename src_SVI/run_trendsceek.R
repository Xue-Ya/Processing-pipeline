library('trendsceek')
library('spatstat')

args = commandArgs()

if (length(args)==0) {
  stop("not enough input", call.=FALSE)
}

count_f <- args[4]
meta_f <- args[5]
out_f <- args[6]

counts <- read.csv(count_f, row.names=1, check.names=F, stringsAsFactors=FALSE)
colData <- read.csv(meta_f, stringsAsFactors=FALSE, row.names=1, check.names=F)


pos2pp_ <- function (pos_mat) 
{
    x.range = range(pos_mat[, 1])
    y.range = range(pos_mat[, 2])
    pp = ppp(x = pos_mat[, 1], y = pos_mat[, 2], x.range, 
        y.range)
    return(pp)
}

gen_null <- function(pp, n.rand = 2000){
    ##Randomize to create null
    ##List of pps with fixed positions but randomized labels

    perm = TRUE #keep the marginal mark dist the same
    pp.perm.list = lapply(1:n.rand, function(j.it, pp, permute){rlabel(pp, permute = permute)}, pp = pp, permute = perm)    
    
    return(pp.perm.list)
}

calc_pp_trendstats <- function(pp, pp.perm.list, alpha_env = 0.05, alpha_nom_early = 0.5, methods = c('Emark', 'markcorr', 'markvario', 'Vmark')){
###Calculate test statistic of obs and null
    
    ##available test stat methods
    fcns = c(Emark, markcorr, markvario, Vmark)
    names(fcns) = c('Emark', 'markcorr', 'markvario', 'Vmark')

    ##subset on selected methods
    fcns = fcns[methods]

    ##calc stats for each selected method
    nsim = length(pp.perm.list)    
    tstat.list = list()
    for(j.fcn.name in names(fcns)){
        j.fcn = fcns[[j.fcn.name]]
        
        print(j.fcn.name)        

        tstat.list[[j.fcn.name]] = calc_pp_trendstats_jmethod(pp, pp.perm.list, j.fcn, alpha_env, alpha_nom_early)
    }

    return(tstat.list)
}

add_subsetstats <- function(local.stats.list, j_localstats, j_nsim_it){
    
    j_nullstats.local.df = localstat2nulldf(j_localstats)            
    
    ##add test-stats from previous subsets
    if(j_nsim_it == 1){
        tstat.local.df = localstat2df(j_localstats)
        nullstats.local.df = j_nullstats.local.df

    }else{
        
        ##add null stats from previous subset
        nullstats.local.df = local.stats.list[['nullstats']]
        nullstats.local.df = cbind(nullstats.local.df, j_nullstats.local.df)

        ##recalculate null-dist based values present in tstat.df, 'r'
        ##and 'obs' are the same regardless of null dist
        tstat.local.df = local.stats.list[['tstat']]
        tstat.local.df[, 'exp'] = apply(nullstats.local.df, 1, mean)
        tstat.local.df[, 'lo.local'] = apply(nullstats.local.df, 1, min)
        tstat.local.df[, 'hi.local'] = apply(nullstats.local.df, 1, max)
    }

    local.stats.list = list(tstat = tstat.local.df, nullstats = nullstats.local.df)

    return(local.stats.list)    
}

localstat2df <- function(tstat){

    ##tstat -> df    
    tstat.df = as.data.frame(tstat)
    cols = colnames(tstat.df)
    cols[which(cols == 'mmean')] = 'exp'
    cols[which(cols == 'lo')] = 'lo.local'
    cols[which(cols == 'hi')] = 'hi.local'
    colnames(tstat.df) = cols

    return(tstat.df)
}


get_envstats <- function(tstat.df, nullstats.df, alpha_env){
    
    ##test-stat exp value can differ between different radii, therefore look at deviation (distance) from the exp and not the test-stat directly
    null.mean = tstat.df[, 'exp']
    null.devs = apply(nullstats.df, 2, function(j.sim, null.mean){abs(j.sim - null.mean)}, null.mean = null.mean)
    ##TBD: variance is scale-dependent, varying with exp, which one may want to account for
    
    null.sorted = t(apply(null.devs, 1, sort)) #sort every r_i
    null.global = sort(apply(null.sorted, 2, max), decreasing = TRUE) #max of every kth ranked value across all r_i

    ##*###
    ##P-value calculation
    ##*###
    nsim = ncol(nullstats.df)
    obs.dev = abs(tstat.df[, 'obs'] - null.mean)
    hi.rank = unlist(lapply(obs.dev, function(j.r, null.global){length(which(null.global >= j.r))}, null.global = null.global))
    hi.rank = hi.rank + 1
    p = hi.rank / (nsim + 1) ##no need to use lo.rank and take *2 since already accounted for two-sided test by taking absolute values

    ##*###
    ##Global envelope
    ##*###
    ##already a two.sided test since took absolute values above
    kth.nullval = floor((nsim + 1) * alpha_env)
    if(kth.nullval == 0){
        kth.nullval = 1
    }
    max.dev = null.global[kth.nullval]

    ##global dev
    hi.global = null.mean + max.dev
    lo.global = null.mean - max.dev

    ##Relative deviation from envelope
    env.rel.dev = obs.dev / max.dev
    
    ##*###
    ##Effect size
    ##*###
    ##relative deviation from null.mean
    mean.rel.dev = obs.dev / null.mean

    ##Store
    envstats = cbind(obs.dev, lo.global, hi.global, env.rel.dev, mean.rel.dev, p)

    return(envstats)
}

localstat2nulldf <- function(tstat){
###get null stats
    
    nullstats.df = as.data.frame(attr(tstat, 'simfuns'))

    ##exclude radii
    nullstats.df.nor = nullstats.df[, setdiff(colnames(nullstats.df), 'r')]

    return(nullstats.df.nor)
}


calc_pp_trendstats_jmethod <- function(pp, pp.perm.list, j.fcn, alpha_env = 0.05, alpha_nom_early = 0.5){
###Test stats for observed and simulated null data

    nsim = length(pp.perm.list)    

    ##alpha_early = alpha * alpha_earlystop_mult. Increase nsim one order of magnitude at a time, e.g. 10, 100, 1000. 
    nsim_order = floor(log10(nsim))
    nsim_vec = 10^(1:nsim_order)
    if(nsim != nsim_vec[length(nsim_vec)]){
        nsim_vec = c(nsim_vec, nsim)
    }
    
    ##loop over increasing nsim
    n_simsplits = length(nsim_vec)
    local.stats.list = list()
    j_nsim_it = 0
    j_p = 0
    alpha_nom_earlystop = 0.5 ## just used for 10 first permutations
    tstat.list = list()
    while(j_p <= alpha_nom_earlystop && j_nsim_it < n_simsplits){
        j_nsim_it = j_nsim_it + 1
        if(j_nsim_it > 1){
            alpha_nom_earlystop = alpha_nom_early
        }
                
        ##select subset from pregenerated null (pp.perm.list)
        ##e.g. 11-100
        if(j_nsim_it == 1){
            jperm_start = 1
        }else{
            jperm_start = nsim_vec[j_nsim_it - 1] + 1
        }            
        jperm_end = nsim_vec[j_nsim_it]
        
        j_pp_perm_list = pp.perm.list[jperm_start:jperm_end]

        
        ##get local (per radius) test stats for subset
        ##This is the time-consuming step
        j_nsim = length(j_pp_perm_list)
        j_localstats = envelope(pp, j.fcn, nsim = j_nsim, simulate = j_pp_perm_list, global = FALSE, nsim2 = 0, savefuns = TRUE, verbose = FALSE)
        
        ##add test-stats from previous subsets
        local.stats.list = add_subsetstats(local.stats.list, j_localstats, j_nsim_it)
        
        ##get point-wise envelope test-stats and p-value
        j_envstats.df = get_envstats(local.stats.list[['tstat']], local.stats.list[['nullstats']], alpha_env) ##tstat.local.df, nullstats.local.df
        
        ##global envelope p-value (min across all radii)
        j_p = min(j_envstats.df[, 'p'])
    }

    if(j_nsim_it < n_simsplits){
        earlystop = 1
    }else{
        earlystop = 0
    }
    
    ##*##    
    ##store stats in list
    ##*##
    
    ##pointwise stats: for every radius value
    tstat.list[['r.stat']] = cbind(local.stats.list[['tstat']], j_envstats.df)
    
    ##global envelope stats
    tstat.list[['max.env.rel.dev']] = max(j_envstats.df[, 'env.rel.dev'])
    tstat.list[['max.rel.dev']] = max(j_envstats.df[, 'mean.rel.dev'])
    tstat.list[['min.pval']] = j_p
    tstat.list[['earlystop']] = earlystop
    tstat.list[['nsim_stop']] = j_nsim_it
    tstat.list[['nsim_max']] = n_simsplits
    
    return(tstat.list)
}

calc_trendstats <- function(jfeat_it, pp, n.rand = 1e4, alpha_env = 0.05, alpha_nom_early = (0.05 * 4) / ifelse(ifelse(is.numeric(pp[['marks']]), length(pp[['marks']]), ncol(pp[['marks']])) >= 500, 10, 1), methods = c('Emark', 'markcorr', 'markvario', 'Vmark')){

    print(jfeat_it)
    
    ##subset on a single mark
    marx = pp[['marks']]
    if(!is.numeric(marx)){
        pp[['marks']] = marx[, jfeat_it]
    }

    ##generate null dist by permuting marks
    ##the same null-dist is used for all methods for comparability
    pp.perm.list = gen_null(pp, n.rand)
    ##TBD: Possibly need to use pool.envelope if running into mem-constraints

    ##calc stats
    tstat.list = calc_pp_trendstats(pp, pp.perm.list, alpha_env, alpha_nom_early, methods)

    return(tstat.list)
}

trendsceek_test <- function(pp, nrand = 1e4, ncores = 1, alpha_env = 0.1 / ifelse(is.numeric(pp[['marks']]), length(pp[['marks']]), ncol(pp[['marks']])), alpha_bh = 0.05, alpha_nom_early = (alpha_bh * 4) / ifelse(ifelse(is.numeric(pp[['marks']]), length(pp[['marks']]), ncol(pp[['marks']])) >= 500, 10, 1)){
    
    ##init parallelization
    bp_param = BiocParallel::MulticoreParam(workers = ncores)

    ##get rstats
    marx = pp[['marks']]
    if(is.numeric(marx)){
        nfeats = 1
        feats = 1
    }else{
        nfeats = ncol(marx)
        feats = colnames(marx)
    }
    tstat_list = BiocParallel::bplapply(1:nfeats, calc_trendstats, BPPARAM = bp_param, pp = pp, n.rand = nrand, alpha_env = alpha_env, alpha_nom_early = alpha_nom_early)
    names(tstat_list) = feats

    ##get supinum stats
    supstats_list = tstat2supstat(tstat_list)
    supstats_wide = supstats_list2wide(supstats_list)
    
    ##store
    trendstat_list = list(tstat = tstat_list, supstats = supstats_list, supstats_wide = supstats_wide)

    ##get sig genes
    sig_list = extract_sig_genes(trendstat_list, alpha_bh)

    trendstat_list[['sig_genes_list']] = sig_list
    
    return(trendstat_list)
}

pp <- pos2pp_(colData)
log.fcn = log10
pp = set_marks(pp, as.matrix(t(counts)), log.fcn = log.fcn)

min.ncells.expr = 3
min.expr = 5
counts_filt = genefilter_exprmat(as.matrix(t(counts)), min.expr, min.ncells.expr)

quantile.cutoff = 0.9 ##filter out the most lowly expressed genes from the fitting
method = 'glm' ##For (robust) linear regression set to 'rlm'
vargenes_stats = calc_varstats(counts_filt, counts_filt, quant.cutoff = quantile.cutoff, method = method)



nrand = 100
ncores = 1

trendstat_list = trendsceek_test(pp, nrand, ncores)

write.csv(trendstat_list[['supstats_wide']], paste0(out_f,"trendsceek.csv"), row.names = TRUE)