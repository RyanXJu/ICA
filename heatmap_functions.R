heatmap(IC_traits_matrix,distfun=function(x) dist(x, method="euclidean"), 
        hclustfun=function(x) hclust(x, method="complete"),
        reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean))


heatmap.2(IC_traits_matrix, col=redgreen(100),
          distfun=function(x) dist(x, method="euclidean"),
          hclustfun=function(x) hclust(x, method="complete"),
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          density.info="none", trace="none", dendrogram=c("row"), 
          symm=F,symkey=T,symbreaks=T, scale="none") 


heatmap(IC_traits_matrix)
heatmap.2(IC_traits_matrix)