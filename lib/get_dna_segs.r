get_dna_segs<-function(nb){
  DL<-list()
  NL<-nb

  for(i in 1:length(nb)) {
    name=nb[[i]]$query.accession
    start<-vector(length=nrow(nb[[i]]))
    end<-vector(length=nrow(nb[[i]]))
    direct<-nb[[i]]$dir
    strand<-ifelse(direct=="fwd",1,-1)

    STRAND<-strand
    STRAND
    if("PF00501" %in% nb[[i]]$PfamID1){
      if(nb[[i]]$dir[nb[[i]]$PfamID1=="PF00501"][1]=="rev"){
        strand<-strand * (-1)
        ord<-sort(1:nrow(nb[[i]]),decreasing=T)

        for(k in 1:length(ord)){
          NL[[i]][k,]<-nb[[i]][(ord[k]),]
          STRAND[k]<-strand[ord[k]]
        }
      }
    }
    
    nb[[i]]<-NL[[i]]
    strand<-STRAND
    pfam<-nb[[i]]$PfamID1
    
    for(j in 1:nrow(nb[[i]])){
      if(j==1){
        start[j]<-nb[[i]]$min[j]
        end[j]<-nb[[i]]$max[j]
      }
      else{
        start[j]<-end[j-1] + 100
        end[j]<-((nb[[i]]$max[j] - nb[[i]]$min[j]) + start[j] + 100)
      }
    }
    
    df1<-data.frame(name=name,
                    start=start,
                    end=end,
                    strand=strand,
                    pfam=pfam)
    dna_seg1<-dna_seg(df1)
    num<-length(dna_seg1$col)
    dna_seg1$col<-"black"
    dna_seg1$fill<-"gray"
    dna_seg1$gene_type<-"arrows"
    DL[[i]]<-dna_seg1
  }
  return(DL)
}


