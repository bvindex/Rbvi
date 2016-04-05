#' Biological Value Index
#' 
#' @description Calculates the Biological Value Index (Sanders, 1960) with modifications from Loya-Salina & Escofet (1990). It also provides a relative index, that allows for comparisons.
#' 
#' @param data an object of class matrix or data.frame. rownames must be samples and colnames must be species (or any other taxa).
#' @param p The cutoff percentage to use, expressed relative to 1 (i.e. for a cutoff of 95% use p = 0.95)
#' 
#' @return results An object of class data.frame containing a column with species name (spp), the scores for each species in each sample, a column (BVI) with the Biological Value Index, and a column (rBVI) with the Relative Biological Value Index.
#' 
#' @seealso bvi_plot
#' 
#' @author Villasenor-Derbez, J.C.
#' 
#' @export

bvi=function(data,p=0.95, other=T){
  
  if (p>1){
    stop("Your cutoff percentage must be between 0 and 1")
  }
  
  if(!is.matrix(data)){
    data=as.matrix(data)
  }
  
  prop=prop.table(data,2) #Calculate matrix of relative abundances by sample
  
  ni=rowSums(data)/sum(data)#Calculate total relative abundances by species
  
  nsp=rep_len(0,dim(data)[2])#vector with zeros of length=number of samples
  
  #The folowing cycle calculates the number of species (nsp) needed to reach the cutoff percentage (p)
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
  #Extract the maximum number of species needed and overwrite nsp
  nsp=max(nsp)
  
  #Pre-demension the matrix of scores so that it is equal to data, but with scores of 0
  scores=data*0
  
  #Cycle iterates across samples (j)
  for (j in seq(1,dim(data)[2])){
    score=unname(prop[,j])        #Extract abundances of sample j
    un=sort(unique(score), decreasing=T)              #Extract unique vector of abundances to correct for ties
    N=nsp                         #N has the maximum number of species needed, and is modified later in the correction for ties
    score_j=score*0               #Set initial score_j = 0 and with size of score
    
    #Cycle iterates across species (i)
    for (i in seq(1,min(length(un)-1,nsp))){
      max_score=max(un)           #max_score contains the maximum abundance to look for during iteration i
      score_j[score==max_score]=N #The score N is assigned to species i (i.e. all those that have an abundance of max?score)
      N=N-sum(score==max_score) #Correction for ties
      un[un==max_score]=0 #Delete the previous values of max_score, so that the next max_score is updatet
    }
    scores[,j]=score_j #Assign scores of sample j to the score matrix
  }
  
  scores=unname(scores)
  BVI=rowSums(unname(scores))
  rBVI=BVI/sum(BVI)*100
  
  results_a=data.frame(rownames(data),
                     scores,
                     BVI,
                     rBVI)
  
  results=results_a[order(-results_a[,dim(scores)[2]+2]),]
  colnames(results)=c("spp",colnames(data),"BVI","%BVI")
  
  if(others){
    results_b=results[1:nsp,]
    results_b[nsp+1,]=colSums(results[nsp+1:dim(results)[1]])
    results=results_b
  }
  
  return(results)
}


