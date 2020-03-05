get_enriched_pfams <- function(noc){

  #Trim the hmm column
  noc <- co
  noct <- noc[, c(1:2,4:(ncol(noc)-4))]

  # Group by the OleC neighborhoods
  nb<-list()
  Cacc<-unique(noct$query.accession)
  for(i in 1:length(Cacc)){
    nb[[i]]<-noct[noct$query.accession==Cacc[i],]
  }
  
  # Look at overall PFAM enrichment
  pfam<-noc$PfamID1
  pfamtab<-data.frame(table(pfam))

  pfamsort<-pfamtab[order(pfamtab$Freq,decreasing=T),]
  pfam10 <- pfamsort[1:12,] 

  # Get PFAM ids
  desc<-unlist(lapply(1:length(pfam10$pfam),function(x){noc$description[grep(pfam10$pfam[x],noc$PfamID1)[1]]}))
  pfam_desc<-data.frame(as.character(pfam10$pfam),pfam10$Freq,desc,stringsAsFactors=F)
  head(pfam_desc)
  pfam_desc$desc<-as.character(pfam_desc$desc)
  pfam_desc$desc[pfam_desc$pfam=="NO PFAM MATCH"]<-"hypothetical protein"

  colnames(pfam_desc)[1:2] <- c("pfam","freq")

  pfam_desc$perc<-(pfam_desc$freq)/length(Cacc)
  pfam_desc$desc[pfam_desc$pfam=="NO PFAM MATCH"]<-"hypothetical protein"
  pfam_desc$desc[pfam_desc$desc==""]<-"hypothetical protein"
  
  
  pfam_desc$pfam<-factor(pfam_desc$pfam,levels=pfam_desc$pfam)
  
  return(pfam_desc)
}