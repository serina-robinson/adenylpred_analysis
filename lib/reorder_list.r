reorder_list<-function(acc,leaves){
  NL<-vector(length=length(acc))
  for(i in 1:length(acc)){
    ind<-grep(names(acc)[i],leaves)
    NL[ind]<-acc[i]
    names(NL[ind])<-leaves[ind]
  }
  return(NL)
}