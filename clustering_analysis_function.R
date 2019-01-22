library(ggplot2)
library(gplots)
library(RColorBrewer)
library(data.table)
library(MGMM)
library(cowplot)
library("rdist")
library(foreach)

# TO DO adapt to dataset different than JASS
compute_silh <- function(dat2)
{
    dat2['Silhouette'] = 0
    Z_score_columns = grep("z_", names(dat2))

    for(i in 1:dim(dat2)[1])
    {
        coli = setdiff(Z_score_columns, which(is.na(dat2[i,Z_score_columns])))
        d = rep(0,length(unique(dat2$cluster)))
        clp = dat2[i, "cluster"]
        unique(dat2$cluster)

        for(cl in unique(dat2$cluster))
        {
            Dc = cdist(dat2[i, coli], dat2[dat2$cluster==cl, coli])
            d[cl] = mean(Dc, na.rm=TRUE)
        }
        mean(Dc, na.rm=TRUE)

        b = min(d[-clp])
        a = d[clp]

        dat2[i,"Silhouette"] = (b - a)/max(a,b)
        dat2[i,"Silhouette"]
    }
    return(dat2)
}

shannon_entropy <- function(x){return(-sum(x*log2(x)))}

add_missing <- function(datf, rm =0.1)
{
    nval = dim(datf)[1] * dim(datf)[2]
    nmissing = floor(rm * nval)

    xi = sample.int(dim(datf)[1], size=nmissing)
    ci = sample.int(dim(datf)[2], size=nmissing, replace=TRUE)

    for(i in 1:nmissing)
    {
        datf[xi[i], ci[i]] = NA
    }
    return(datf)
}


align_zscores <- function(dat2)
{
    # the variable with the most variance will be set
    # to be always positive and other variable are flip
    # according to it

    reference_variable = which.max(sapply(dat2[, grepl("z_", names(dat2))], var,na.rm=TRUE))
    ids_d = which(dat2[,reference_variable] < 0)
    print( paste("SNP aligned on", names(dat2)[reference_variable]))
    for(c in names(dat2)){
        dat2[ids_d, c] = -dat2[ids_d, c]
    }

    return(dat2)
}

initialize_at_random <- function(dat2, cn)
{
    M0 = list()
    S0 = list()
    # If there is too few complete sample do a randomize centroid
    # initialization
    ids = sample(dim(dat2)[1],cn, replace=FALSE)
    med1 = apply(dat2, 2, median, na.rm=TRUE) # median input
    for(i in 1:cn)
    {
        M0[[i]] = as.numeric(dat2[ids[i],])
        S0[[i]] = diag(dim(dat2)[2])
        idna = which(is.na(M0[[i]]))
        M0[[i]][idna] = med1[idna]
        # choose a random data point
    }
    pi0 = rep(1.0/cn, cn)
    return(list(M0=M0, S0=S0, pi0=pi0))
}


#dat2 must be a dataframe.....
compute_data_cluster <- function(dat2,cn)
{
    print(paste("Compute", cn, "clusters" ))
    dat2 = dat2[, grepl( "z_", names(dat2))]
    # Apply clustering
    print(sum(complete.cases(dat2)))
    print(cn*7)
    if(sum(complete.cases(dat2)) < (dim(dat2)[2]*cn*5))
    {
        print("initialize at random")
        L = initialize_at_random(dat2, cn)
        Mt = fit.MNMix(as.matrix(dat2), M0=L$M0, S0=L$S0,
         L$pi0, k=cn, parallel=FALSE, report=TRUE, maxit=500)
    }
    else{
        print("initialize with k-means")
        Mt = fit.MNMix(as.matrix(dat2), k=cn, parallel=FALSE, report=TRUE, maxit=500)
    }
    print(Mt)

    dat2["cluster"]= Mt@Assignments$Assignment
    dat2["entropy"] = apply(Mt@Responsibilities[,-1], 1, shannon_entropy)
    dat2 = compute_silh(dat2)
    return(list( clustered_data = dat2, MNMmix=Mt))
}

compute_perf_by_k <- function(dat2,
    nclust=2:10,
    save_clustering=FALSE, save_model=FALSE,
     prefix_file='./Cluster')
{
    n = dim(dat2)[1]
    m = dim(dat2)[2]

    data_clust_perf = foreach( i = nclust, .combine="rbind") %do% {
        tryCatch({

        cl_res = compute_data_cluster(dat2, i)
        #dat2 = cl_res$clustered_data
        Mt = cl_res$MNMmix

        if(save_clustering){write.csv(cl_res$clustered_data, file=paste0(prefix_file, i, ".csv"))}
        if(save_model){saveRDS(Mt, file=paste0(prefix_file, i, "_MNM_model.rds"))}

        result = c(k=i, Objective = Mt@Objective,
            Silhouette =mean(cl_res$clustered_data$Silhouette),
            BIC = (i*m*log(n) - 2*Mt@Objective))
            return(result)
        }, error = function(e){
            print(paste('Error', e, "Occured for task", i,
            "\n returning a NULL vector"))
            return (rep(NULL, 4))
        }
      )
    }
    data_clust_perf = as.data.frame(data_clust_perf)
    names(data_clust_perf) = c("k", 'Objective', 'Silhouette', "BIC")

    return(data_clust_perf)
}

### Plotting function
compute_perf_by_k_bootstrap <- function(dat2,
    maxclust=10, repetition=10)
{

    n = dim(dat2)[1]
    m = dim(dat2)[2]
    nrep = ((10-2)*repetition)

    data_clust_perf = foreach( i = 1:nrep, .combine="rbind") %do% {
        tryCatch({

        nc = i %% (maxclust-2) + 2

        print(nc)
        print(paste("n cluster : ", nc))

        ids = sample(n, round(n*0.8), replace=TRUE)
        print(head(dat2))
        cl_res = compute_data_cluster(dat2[ids,], nc)
        Mt = cl_res$MNMmix

        result = c(id_rep = i, k=nc, Objective = Mt@Objective,
            Silhouette =mean(cl_res$clustered_data$Silhouette),
            BIC = (nc*m*log(n) - 2*Mt@Objective))

            return(result)
        }, error = function(e){
            print(paste('Error', e, "Occured for task", i,
            "\n returning a NULL vector"))
            return (rep(NULL, 4))
        }
    )
    }

    data_clust_perf = as.data.frame(data_clust_perf)
    names(data_clust_perf) = c("Id", "k", 'Objective', 'Silhouette', "BIC")

    return(data_clust_perf)
}


 optimal_k_plot <- function(data_clust_perf)
 {
     perf1 = qplot(k, Objective, data=data_clust_perf, colour=I('royalblue'), geom=c("point", "line"), size=I(1.5))
     k_best = which.min(data_clust_perf$Objective) + 1
     Objective_best = min(data_clust_perf$Objective)
     perf1 = perf1 + annotate("point",x=I(k_best), y=I(Objective_best), color=I('firebrick1'), size=5, shape=17)

     perf2 = qplot(k, BIC, data=data_clust_perf, colour=I('royalblue'), geom=c("point", "line"), size=I(1.5))
     k_best = which.max(data_clust_perf$BIC) + 1
     Objective_best = max(data_clust_perf$BIC)
     perf2 = perf2 + annotate("point",x=I(k_best), y=I(Objective_best), color=I('firebrick1'), size=5, shape=17)

     perf3 = qplot(k, Silhouette, data=data_clust_perf, colour=I('royalblue'),
     geom=c("point", "line"), size=I(1.5))
     k_best = which.max(data_clust_perf$Silhouette) + 1
     Objective_best = max(data_clust_perf$Silhouette)
     perf3 = perf3 + annotate("point",x=I(k_best), y=I(Objective_best), color=I('firebrick1'), size=5, shape=17)

     return(plot_grid(perf1, perf2, perf3, nrow=1, labels=c("a","b","c")))
 }

#
optimal_k_plot_bootstrap <- function(data_clust_perf)
{
    dat_clust = as.data.table(data_clust_perf)

    mean_perf = dat_clust[,lapply(.SD, mean), by=k]
    sd_perf = dat_clust[,lapply(.SD, sd), by=k]

    setkey(mean_perf, "k")
    setkey(sd_perf, "k")

    dat_perf_total = merge(mean_perf, sd_perf, suffixes=c(".mean", ".sd"))

    p1 <- ggplot(dat_perf_total, aes(x=k, y=Objective.mean)) +
      geom_line(color="royalblue") +
      geom_point(color="royalblue")+  geom_jitter(data = dat_clust, width=0.25, alpha=0.5, color="darkgray", aes(x=k, y=Objective)) +
      geom_errorbar(aes(ymin=Objective.mean-Objective.sd, ymax=Objective.mean+Objective.sd), width=0.2,
                     position=position_dodge(0.05), color="royalblue") +xlab('number of cluster')

    p2<- ggplot(dat_perf_total, aes(x=k, y=BIC.mean)) +
    geom_line(color="royalblue") +
    geom_point(color="royalblue")+  geom_jitter(data = dat_clust, width=0.25, alpha=0.5, color="darkgray", aes(x=k, y=BIC))+
    geom_errorbar(aes(ymin=BIC.mean-BIC.sd, ymax=BIC.mean+BIC.sd), width=0.2,
              position=position_dodge(0.05), color="royalblue") +xlab('number of cluster')

    ymin = min(0.1, min(dat_clust$Silhouette))
    ymax = max(0.3, max(dat_clust$Silhouette))

      p3<- ggplot(dat_perf_total, aes(x=k, y=Silhouette.mean)) +
        geom_line(color="royalblue") + ylim(c(ymin, ymax)) +
        geom_point(color="royalblue")+ geom_jitter(data = dat_clust, width=0.25, alpha=0.5, color="darkgray", aes(x=k, y=Silhouette))+ geom_errorbar(aes(ymin=Silhouette.mean-Silhouette.sd, ymax=Silhouette.mean+Silhouette.sd), width=0.2, position=position_dodge(0.05), color="royalblue") +xlab('number of cluster')


     return(plot_grid(p1, p2, p3, nrow=1, labels=c("a","b","c")))
}
