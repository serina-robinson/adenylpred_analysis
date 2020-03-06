make_cluster <- function(shrt, thresh, bynum, namvec, pal) { # input is an all vs. all BLAST table 

  # First remove all comparisons between identical proteins (closed loop nodes)
  noprs <- shrt[!shrt$prot1==shrt$prot2,]
  noprs <- noprs[order(noprs$eval),]
  
  source("lib/color_clustering_diagram.r")
  meta <- color_nodes(noprs, namvec, pal)
  
  # Make a unique data frame 
  metu<-unique(meta)

  # Set a similarity threshold 
  net <- noprs[noprs$eval<(as.numeric(paste0("1.00e-",thresh))),]
  
  metu2 <- metu[order(metu$size),]

  # Simplify the graph
  g <- simplify(graph.data.frame(net, vertices=metu2, directed = FALSE))
  g <- delete.vertices((g), degree(g) < 1) # remove singletons
  
  # Append metadata
  V(g)$color <- V(g)$color
  V(g)$size <- V(g)$size
  
  l <- layout_components(g)
  
  # Make a legend
  colors<-unique(metu$color)
  proteins<-unique(metu$fam)
  clegend <- data.frame(colors, proteins)
  
  return(list(net=net, g=g, legend = clegend))
}
