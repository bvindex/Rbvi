#' 
#' 
#' 
#' 
#' 

bvi=function(data,p=0.95, r=FALSE){
  
  if(!is.matrix(data)){
    data=as.matrix(data)
  }
  
  prop=prop.table(data,2) #Calculate matrix of relative abundances by sample
  
  ni=rowSums(data)/sum(data)#Calculate total relative abundances by species
  
  nsp=rep_len(0,dim(data)[2])#vector with zeros of length=number of samples
  
  for (j in seq(1,dim(data)[2])){
    i=1
    nij=sort(prop[,j],decreasing=T)
    ni_p=nij[i]
    
    while(ni_p<p){
      ni_p=ni_p+nij[i+1]
      i=i+1
    }
    nsp[j]=i
  }
  nsp=max(nsp)
  
  scores=matrix(0,nrow=dim(data)[1], ncol=dim(data)[2])
  
  for (j in seq(1,dim(data)[2])){
    score=unname(prop[,j])
    un=sort(unique(score))
    N=nsp
    score_j=score*0
    
    for (i in seq(1,min(length(un)-1,nsp))){
      max_score=max(un)
      score_j[score==max_score]=N
      N=N-sum(score==max_score)
      un[un==max_score]=0
    }
    scores[,j]=score_j
  }
  
  scores=sort(rowSums(scores), decreasing=T)
  others=sum(scores[(nsp+1):length(scores)])
  BVI=c(scores[1:nsp],others)

  results=data.frame(spp=c(rownames(data[1:nsp,]),"Others"),
                     BVI)
  
    
  if (r){
    results[,3]=results$BVI/sum(results$BVI)*100
    colnames(results)=c("spp","BVI", "%BVI")
  }
  
  return(results)
}


