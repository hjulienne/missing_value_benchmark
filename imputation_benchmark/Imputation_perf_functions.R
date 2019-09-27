library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)
### Plotting routines ####
lg_style = theme(legend.position="top",
            legend.title = element_text(colour="black", size=10,
                                face="bold"))

impute_by_fun <- function(X, imp_fun =mean )
{

  for(var in names(X))
  {
    X[is.na(X[var]), var] = imp_fun(X[,var], na.rm=TRUE)
  }
  return(X)
}

plot_no_cluster <- function(X_ref, X_test)
{
    p1 = ggplot(X_ref, aes(x=y1, y=y2)) + geom_point() + geom_density2d()
    p2 = ggplot(X_ref, aes(x=y1, y=y2)) + geom_point() + geom_density2d()
    plot_grid(p1,p2)
}

plot_against_reference <- function(X_ref, X_test)
{
    p1 = ggplot(X_ref, aes(x=y1, y=y2,color=component)) + geom_point() + geom_density2d()
    p1 = p1 + lg_style + guides(color=guide_legend(title="component"))
    p2 = ggplot(X_ref, aes(x=y1, y=y2)) + geom_point(aes(color= as.factor(X_test$cluster))) + geom_density2d()
    p2 = p2 + lg_style + guides(color=guide_legend(title="cluster"))

    plot_grid(p1,p2)
}

plot_reconstruction_against_reference <- function(X_ref, X_test)
{
    p1 = ggplot(X_ref, aes(x=y1, y=y2,color=component)) + geom_point() + geom_density2d()
    p1 = p1 + lg_style + guides(color=guide_legend(title="component"))
    p2 = ggplot(X_test, aes(x=y1, y=y2)) + geom_point(aes(color= as.factor(X_test$cluster))) + geom_density2d()
    p2 = p2 + lg_style + guides(color=guide_legend(title="cluster"))

    plot_grid(p1,p2)
}

plot_performance_dashboard <- function(dat_perf)
{
  dat_long= melt(dat_perf, id.vars= "X")

  p = ggplot(dat_long, aes(x=X, y = value, color=variable)) + geom_line(lwd=1.25) + geom_point(size=2.5)
  p = p + xlab("missing data ratio") + ylab("Adjusted Rand Index") + theme_linedraw()
  p = p + lg_style + scale_colour_brewer(palette = "Set1")# brewer.pal(6,"Set1")
  p = p + guides(color=guide_legend(title="method"))
  return(p)
}

plot_performance_adj_error <- function(folder, pattern)
{
  dat_list = list()
  file_list = list.files(folder, pattern= pattern)
  file_list
  for(d in file_list){
   dat_list[[d]] = read.csv(paste0(folder,d))
  }

  N = length(file_list)
  data_bm = bind_rows(dat_list)
  data_bm = as.data.table(data_bm)

  DM = data_bm[, lapply(.SD, mean), by=X]
  Dsd = data_bm[, lapply(.SD, sd), by=X]
  did_compute = function(x){sum(x >0)}
  N_simu = data_bm[, lapply(.SD, did_compute),by=X]

  dat_long= as.data.frame(melt(DM, id.vars= "X"))
  dat_long_error = as.data.frame(melt(Dsd, id.vars= "X"))
  dat_long_count = as.data.frame(melt(N_simu, id.vars="X"))

  is_computation_performed = which(dat_long$value > 0 & dat_long_count$value > 1)
  dat_long = dat_long[is_computation_performed,]
  dat_long_error = dat_long_error[is_computation_performed,]
  dat_long_count = dat_long_count[is_computation_performed,]

  dat_long_error["ymin"] = dat_long["value"] - dat_long_error["value"] / (dat_long_count["value"]^0.5)
  dat_long_error["ymax"] = dat_long["value"] + dat_long_error["value"] / (dat_long_count["value"]^0.5)

  p = ggplot(dat_long, aes(x=X, y = value, color=variable)) + geom_line(lwd=1.25) + geom_point(size=2.5)
  p = p + xlab("missing data ratio") + ylab("Adjusted Rand Index") + theme_linedraw()
  p = p  + lg_style+ scale_colour_brewer(palette = "Set1")# brewer.pal(6,"Set1")+ lg_style
  p = p + guides(color=guide_legend(title="method"))
  p = p + geom_errorbar(data=dat_long_error, mapping=aes(x=X, ymin=ymin, ymax=ymax, color=variable))

  return(p)
  #
  # file_list = list.files(folder, pattern= pattern)
  # for(d in file_list)
  # {
  #  dat_list[[d]] = read.csv(paste0(folder,d))
  # }
  # N = length(file_list)
  # data_bm = bind_rows(dat_list)
  # data_bm = as.data.table(data_bm)
  #
  # DM = data_bm[, lapply(.SD, mean), by=X]
  # Dsd = data_bm[, lapply(.SD, sd), by=X]
  #
  # dat_long= as.data.frame(melt(DM, id.vars= "X"))
  # dat_long_error = as.data.frame(melt(Dsd, id.vars= "X"))
  # dat_long_error["ymin"] = dat_long["value"] - dat_long_error["value"] / (N^0.5)
  # dat_long_error["ymax"] = dat_long["value"] + dat_long_error["value"] / (N^0.5)
  #
  # p = ggplot(dat_long, aes(x=X, y = value, color=variable)) + geom_line(lwd=1.25) + geom_point(size=2.5)
  # p = p + xlab("missing data ratio") + ylab("Adjusted Rand Index") + theme_linedraw()
  # p = p + lg_style + scale_colour_brewer(palette = "Set1")# brewer.pal(6,"Set1")
  # p = p + guides(color=guide_legend(title="method"))
  # p = p + geom_errorbar(data=dat_long_error, mapping=aes(x=X, ymin=ymin, ymax=ymax, color=variable))
  # return(p)
}



# Create function
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
