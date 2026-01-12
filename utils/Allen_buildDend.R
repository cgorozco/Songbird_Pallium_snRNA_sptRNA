
###---- ALLEN INSTITUTE CLUSTERGRAM
build_dend <- function(cl.dat, cl.cor=NULL, l.rank=NULL, l.color=NULL, nboot=100, ncores=1)
{
  require(dendextend)
  require(dplyr)
  if(is.null(cl.cor)){
    cl.cor = cor(cl.dat)
  }
  pvclust.result=NULL
  if(nboot > 0){
    require(pvclust)
    parallel= FALSE
    if(ncores > 1){
      parallel = as.integer(ncores)
    }
    pvclust.result <- pvclust::pvclust(cl.dat, method.dist = "cor" ,method.hclust = "average", nboot=nboot, parallel=parallel)
    dend = as.dendrogram(pvclust.result$hclust)
    dend = label_dend(dend)$dend
    dend = dend %>% pvclust_show_signif_gradient(pvclust.result, signif_type = "bp", signif_col_fun=colorRampPalette(c("white","gray","darkred","black")))
    #%>% pvclust_show_signif(pvclust.result, signif_type="bp", signif_value=c(2,1))
  }
  else{
    cl.hc = hclust(as.dist(1-cl.cor),method="average")      
    dend = as.dendrogram(cl.hc)
  }
  dend = dend %>% set("labels_cex", 0.7)
  if(!is.null(l.color)){
    dend = dend %>% set("labels_col", l.color[labels(dend)])
  }
  dend = dend %>% set("leaves_pch", 19) %>% set("leaves_cex", 0.5)
  if(!is.null(l.color)){
    dend = dend %>% set("leaves_col", l.color[labels(dend)])
  }
  if(!is.null(l.rank)){
    dend =reorder_dend(dend,l.rank)
  }
  return(list(dend=dend, cl.cor=cl.cor, pvclust.result=pvclust.result))
}



label_dend <- function(dend,n=1)
{  
  if(is.null(attr(dend,"label"))){
    attr(dend, "label") =paste0("n",n)
    n= n +1
  }
  if(length(dend)>1){
    for(i in 1:length(dend)){
      tmp = label_dend(dend[[i]], n)
      dend[[i]] = tmp[[1]]
      n = tmp[[2]]
    }
  }
  return(list(dend=dend, n))
}


