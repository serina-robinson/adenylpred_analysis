color_segs <- function(raw_segs, pfams){

  pal <- c("gray", "#FFB341", "#92D050", "#F37FCC", "dodgerblue", "#E31A1C", "#CAB2D6",
           "#6A3D9A", "#FFFF99", "#B15928", "#34A02C", "black")

  for(i in 1:length(raw_segs)){
    for(j in 1:length(pfams)){
      ind<-grep(pfams[j],x = raw_segs[[i]]$pfam)
      raw_segs[[i]]$col[ind]<-"black"
      raw_segs[[i]]$fill[ind]<-pal[j]
    }
  }  
  return(raw_segs)
}

