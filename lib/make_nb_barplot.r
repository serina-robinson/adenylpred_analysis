make_nb_barplot<-function(noc){

  #Trim the hmm column
  noc<-co
  noct<-noc[,c(1:2,4:(ncol(noc)-4))]

  #Group by the OleC neighborhoods
  nb<-list()
  Cacc<-unique(noct$query.accession)
  for(i in 1:length(Cacc)){
    nb[[i]]<-noct[noct$query.accession==Cacc[i],]
  }
  
  #Look at overall PFAM enrichment
  pfam<-noc$PfamID1
  pfamtab<-data.frame(table(pfam))

  pfamsort<-pfamtab[order(pfamtab$Freq,decreasing=T),]
  pfam10<-pfamsort[1:12,] 

  #Get PFAM ids
  desc<-unlist(lapply(1:length(pfam10$pfam),function(x){noc$description[grep(pfam10$pfam[x],noc$PfamID1)[1]]}))
  pfam_desc<-data.frame(as.character(pfam10$pfam),pfam10$Freq,desc,stringsAsFactors=F)
  head(pfam_desc)
  pfam_desc$desc<-as.character(pfam_desc$desc)
  pfam_desc$desc[pfam_desc$pfam=="NO PFAM MATCH"]<-"hypothetical protein"

  colnames(pfam_desc)[1:2]<-c("pfam","freq")
  #return(pfam_desc)

  pfam_desc$perc<-(pfam_desc$freq)/length(Cacc)
  pfam_desc$desc[pfam_desc$pfam=="NO PFAM MATCH"]<-"hypothetical protein"
  pfam_desc$desc[pfam_desc$desc==""]<-"hypothetical protein"
  
  
  pfam_desc$pfam<-factor(pfam_desc$pfam,levels=pfam_desc$pfam)
  source("src/theme.publication.r")
  #pdf("output/cysteine_hydrolase_neighborhoods.pdf",width=25,height=15)
  #palette(colorRampPalette(colors=brewer.pal(12,"Paired"))(12))
  #pal<-palette(colorRampPalette(colors=brewer.pal(12,"Paired"))(12))
  pal<-c("gray","#FB9A99","#B2DF8A","dodgerblue","#33A02C","#E31A1C","#FDBF6F",
         "darkorange1","#CAB2D6","#6A3D9A","#FFFF99","#B15928")
  palette(c("gray","#FB9A99","#B2DF8A","dodgerblue","#33A02C","#E31A1C","#FDBF6F",
            "darkorange1","#CAB2D6","#6A3D9A","#FFFF99","#B15928"))
  gr<-ggplot(data=pfam_desc,aes(y=perc,x=pfam))+
    geom_bar(stat="identity",fill=pal) +
    theme.publication(pfam) +
    scale_x_discrete(name ="PFAM ID") +
    scale_y_continuous(name = "Gene frequency") +
    geom_text(aes(label=round(perc,2)),size=12,vjust=-0.25)
  #dev.off()
  
  pdf(paste0("output/",length(Cacc),"_RODEO_barplot.pdf"),width=25,height=15)
  ggplot(data=pfam_desc,aes(y=perc,x=pfam))+
    geom_bar(stat="identity",fill=pal) +
    theme.publication(pfam) +
    scale_x_discrete(name ="PFAM ID") +
    scale_y_continuous(name = "Gene frequency") +
    geom_text(aes(label=round(perc,2)),size=12,vjust=-0.25)
  dev.off()
  
  pdf(paste0("output/brasiliensis_RODEO_legend.pdf"),width=10,height=10)
  plot.new()
  palette(colorRampPalette(colors=brewer.pal(12,"Paired"))(12))
  l<-legend("bottom",legend=pfam_desc$desc,fill=pal,
         border=FALSE, bty="n", title = "Enzyme")
  l
  dev.off()

  return(list(gr,pfam_desc,l))
}