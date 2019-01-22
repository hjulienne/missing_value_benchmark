#
# two JASS analysis with two trait to enable vizualisation
#
library(doParallel)
library(RColorBrewer)
library(UpSetR)
library(naniar)

source('./clustering_analysis_function.R')
jass_group = c("Composite_complete")#"Immunity_complete_bivariate",
                #"Immunity_complete", "Kidney_faillure_complete")

Nmissing <- function(x){return(sum(is.na(x)))}
setwd("/mnt/atlas/GGS/__PROJECT_GMM_with_missing")
for(id_tr in jass_group)
{
    print("Launching clustering analysis for")
    print(id_tr)

    tryCatch(
    {
        fi = paste0("./datasets/Zscores/", id_tr,".csv")

        Z_data = fread(fi, index="snp_ids")#index='SNP')

        print(Z_data[,grepl("z_", names(Z_data)), with=FALSE][,lapply(.SD, var, na.rm=TRUE)])
        Z_data[,grepl("z_", names(Z_data)), with=FALSE][,lapply(.SD, Nmissing)]


        dat2 = as.data.frame(Z_data[,grepl("z_", names(Z_data)), with=FALSE])
        row.names(dat2) = Z_data$SNP

        #too_much_missing <- which(apply(is.na(dat2), 2, mean) > 0.8)
        # Apply clustering

        #dat2 = dat2[apply(is.na(dat2), 1,mean) < 0.5, ]
        head(dat2)
        considered_traits = names(dat2)
        considered_snps = row.names(dat2)

      #  sqrt_sample_size = sqrt(GWAS_labels[names(dat2), "Nsample"])
      #  sqrt_sample_size

      #  norm_data <- function(x){
      #      print(x)
      #      return(x/sqrt_sample_size)
      #  }
        #dat2 = apply(dat2, 1, norm_data)

        #dat2 = as.data.frame(t(dat2))

        print(paste("considered_trait:", considered_traits))
        fraction_missing = 1-mean(is.na(dat2)) #sum(Z_data[,grepl("z_", names(Z_data)), with=FALSE][,lapply(.SD, Nmissing)]) / (dim(dat2)[1]*dim(dat2)[2])

        print(paste("Fraction of complete cases :", sum(complete.cases(dat2)) / dim(dat2)[1] ))
        print(paste("Fraction of missing observation :", fraction_missing))

        dat2 = align_zscores(dat2)
        print(dim(dat2))
        #dat2 = scale(dat2)

        row.names(dat2) = considered_snps
        names(dat2) = considered_traits
        dat2 = as.data.frame(dat2)

        registerDoParallel(cores=10)
        head(dat2)
        dat2[mapply(is.infinite, dat2)] <- NA
        perf_data = compute_perf_by_k_bootstrap(dat2,
             maxclust=9,
            repetition=100)

        write.csv(perf_data, paste0("./clustering_performances/", id_tr, "_Zscore.csv"))
        png(paste0("./criteria_plots/optimal_k_plot_", id_tr, "_Zscore.png" ),
         width=1200, height=400)
            print(optimal_k_plot_bootstrap(perf_data))
        dev.off()

    },
    error = function(e){
        print(paste('Error', e, "Occured for task", id_tr))
    }
    )
}
