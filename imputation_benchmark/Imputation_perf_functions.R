




library(ggplot2)
library(cowplot)



lg_style = theme(legend.position="top",
            legend.title = element_text(colour="black", size=10,
                                face="bold"))



plot_against_reference <- function(X_ref, X_test)
{
    p1 = ggplot(X_ref, aes(x=y1, y=y2,color=component)) + geom_point() + geom_density2d()
    p1 = p1 + lg_style + guides(color=guide_legend(title="component"))
    p2 = ggplot(X_ref, aes(x=y1, y=y2)) + geom_point(aes(color= as.factor(X_test$cluster))) + geom_density2d()
    p2 = p2 + lg_style + guides(color=guide_legend(title="cluster"))

    return(plot_grid(p1,p2))
}


# Create function
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
