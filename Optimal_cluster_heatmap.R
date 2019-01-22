#
# two JASS analysis with two trait to enable vizualisation
#
library(doParallel)
library(RColorBrewer)
library(UpSetR)
library(naniar)
library(colorspace)
source('./clustering_analysis_function.R')

jass_group = c("Composite_complete")#,"Immunity_complete_bivariate",
                #"Immunity_complete", "Kidney_faillure_complete")
optimal_k_group = c(6,2,5,5)

Nmissing <- function(x){ return(sum(is.na(x)))}
id_tr <- jass_group[1]

for(id_i in 1:1)
{
    print("Launching clustering analysis for")
    print(id_tr)

    id_tr = jass_group[id_i]
    optimal_k = optimal_k_group[id_i]

    tryCatch(
    {
        fi = paste0("./datasets/Zscores/", id_tr,".csv")

        Z_data = fread(fi, index='snp_ids')

        print(Z_data[,grepl("z_", names(Z_data)), with=FALSE][,lapply(.SD, var, na.rm=TRUE)])
        Z_data[,grepl("z_", names(Z_data)), with=FALSE][,lapply(.SD, Nmissing)]
        dim(Z_data)

        dat2 = as.data.frame(Z_data[,grepl("z_", names(Z_data)), with=FALSE])
        row.names(dat2) = Z_data$SNP

        dat2 = dat2[apply(is.na(dat2), 1,mean) < 0.5, ]
        considered_traits = names(dat2)
        considered_snps = row.names(dat2)

        print(paste("considered_trait:", considered_traits))
        fraction_missing = 1-mean(is.na(dat2))

        print(paste("Fraction of complete cases :", sum(complete.cases(dat2)) / dim(dat2)[1] ))
        print(paste("Fraction of missing observation :", fraction_missing))

        dat2 = align_zscores(dat2)
        dat2[mapply(is.infinite, dat2)] <- NA

        getwd()
        row.names(dat2) = considered_snps
        names(dat2) = considered_traits

        dat2 = as.data.frame(dat2)

        sum(!complete.cases(dat2[, grepl("z_", names(dat2))]))
        dat2 = compute_data_cluster(dat2[, grepl("z_", names(dat2))], optimal_k)
        dat2 = dat2$clustered_data
        entropy_treshold = 0.75

        class(dat2)
        pdf(paste0("./heatmaps/entropy_distribution_", id_tr,".pdf"))
            p = ggplot(dat2, aes(x=entropy)) + geom_histogram(fill="royalblue")
            print(p +geom_vline(aes(xintercept=entropy_treshold)))
        dev.off()

        percent_well_clustered = mean(dat2$entropy < entropy_treshold)*100
        print(paste("Filtering Snp with high entropy:", round(percent_well_clustered,2) ,"of the dataset retained"))

        dat2 = dat2[dat2$entropy < entropy_treshold,]
        dat2 = dat2[!is.na(dat2[,dim(dat2)[2]]),]
        bks = c(-30,-10,-5, -3, -1, 0, 1,  3, 5,10, 30)
        col_heat = diverge_hsv(10)

        write.csv(dat2, paste0("./Final_clusters/", id_tr, "_k_", optimal_k,".csv"))

        head(dat2[order(dat2$cluster), grepl("z_", names(dat2))])

        pdf(paste0("./heatmaps/heatmap_", id_tr,"_Z_score.pdf" ), width=8, height=24)
            heatmap.2(as.matrix(dat2[order(dat2$cluster), grepl("z_", names(dat2))]),
            breaks = bks,
            margins=c(10,10),
            trace="none",
            na.color = "white", scale="none",
            col = col_heat, Rowv=FALSE, Colv=FALSE,
            RowSideColors = brewer.pal(12, "Set3")[sort(dat2$cluster)],
            cexCol=0.8)
        dev.off()


        dat2 = as.data.table(dat2)
        cols = names(dat2)[grepl("z_", names(dat2))]
        dat2[, lapply(.SD, mean, na.rm=TRUE), .SDcols=cols, by=cluster]

        dat_centroid = dat2[, lapply(.SD, mean, na.rm=TRUE), .SDcols=cols, by=cluster]
        dat_centroid = dat_centroid[dat_centroid$cluster, ]

        write.csv(dat_centroid, paste0("./Final_clusters/", id_tr, "_centroid_k_", optimal_k,".csv"))

        pdf(paste0("./heatmaps/heatmap_centroid_", id_tr,"_Z_score.pdf" ), width=8, height=6)
        heatmap.2(as.matrix(dat_centroid[, grepl("z_", names(dat_centroid)), with=FALSE]),
                breaks = bks,
                margins=c(10,10),
                trace="none",
                na.color = "white", scale="none",
                col = col_heat, Rowv=FALSE, Colv=FALSE,
                RowSideColors = brewer.pal(12, "Set3")[dat_centroid$cluster],
              cexCol=0.8)
        dev.off()
    },
    error = function(e){
        print(paste('Error', e, "Occured for task", id_tr))
    }
    )
}
