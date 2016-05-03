#' Biological Value Index
#' 
#' @description Calculates the Biological Value Index (Sanders, 1960) with modifications from Loya-Salina & Escofet (1990). It also provides a relative index, that allows for comparisons.
#' 
#' @param data an object of class matrix or data.frame. rownames must be samples and colnames must be species (or any other taxa).
#' @param p The cutoff percentage to use, expressed relative to 1 (i.e. for a cutoff of 95 percent use p = 0.95)
#' @param other a logical that indicates if species that do not contribute to p should be added together ina category "other". If FALSE, they are displayed individually.
#' 
#' @return results An object of class data.frame containing a column with species name (spp), the scores for each species in each sample, a column (BVI) with the Biological Value Index, and a column (rBVI) with the Relative Biological Value Index. The results can be directly passed to bvi_plot to get a graphical representation.
#' 
#' @seealso bvi_plot
#' 
#' @author Villasenor-Derbez, J.C.
#' 
#' @export

bvi=function(data,p=0.95, other=T){
  
  #Test that the cutoff percentage is between 0 and 1
  if (p>1){
    stop("Your cutoff percentage must be between 0 and 1")
  }
  
  #Test that the data passed is a matrix, else it transforms it
  if(!is.matrix(data)){
    data=as.matrix(data)
  }
  
  prop=prop.table(data,2)       #Calculate matrix of relative abundances by sample
  
  ni=rowSums(data)/sum(data)    #Calculate total relative abundances by species
  
  nsp=rep_len(0,dim(data)[2])   #vector with zeros of length=number of samples
  
  #The folowing cycle calculates the number of species (nsp) needed to reach the cutoff percentage (p)
  for (j in seq(1,dim(data)[2])){
    i=1                             #Re-set i at 1
    nij=sort(prop[,j],decreasing=T) #Sort the column j in a decreasing order
    ni_p=nij[i]                     #ni_p will be the cumulative relative abundance. Here it is the largest value of the sorted column
    
    #This conditional keeps adding species' relative abundance to ni_p until the cutoff percentage p is reached
    while(ni_p<p){
      ni_p=ni_p+nij[i+1]            #ni_p increases with the addition of every species [i+1]
      i=i+1                         #i keeps track of the number of species that have been added to reach p
    }
    nsp[j]=i                   #For sample j, the number of species to reach p is i
  }
  
  #Extract the maximum number of species needed and overwrite nsp
  nsp=max(nsp)
  
  #Pre-demension the matrix of scores so that it is equal to data, but with scores of 0
  scores=data*0
  
  #Cycle iterates across samples (j)
  for (j in seq(1,dim(data)[2])){
    score=unname(prop[,j])                 #Extract abundances of sample j
    un=sort(unique(score), decreasing=T)   #Extract unique vector of abundances to correct for ties
    N=nsp                                  #N has the maximum number of species needed, and is modified later in the correction for ties
    score_j=score*0                        #Set initial score_j = 0 and with size of score
    
    #Cycle iterates across species (i)
    for (i in seq(1,min(length(un)-1,nsp))){
      max_score=max(un)           #max_score contains the maximum abundance to look for during iteration i
      score_j[score==max_score]=N #The score N is assigned to species i (i.e. all those that have an abundance of max?score)
      N=N-sum(score==max_score)   #Correction for ties
      un[un==max_score]=0         #Delete the previous values of max_score, so that the next max_score is updatet
    }
    scores[,j]=score_j            #Assign scores of sample j to the score matrix
  }
  
  scores=unname(scores)           #Unname scores so that BVI does not carry with the names
  BVI=rowSums(unname(scores))     #Calculates BVI by adding scores for each species across samples
  rBVI=BVI/sum(BVI)*100
  
  #Build preeliminar results into a data.frame
  results_a=data.frame(rownames(data),
                     scores,
                     BVI,
                     rBVI,
                     stringsAsFactors=FALSE)
  
  #Results finalized
  results=results_a[order(-results_a[,dim(scores)[2]+2]),]
  colnames(results)=c("spp",colnames(data),"BVI","%BVI")
  
  #If other=TRUE, adds species that do not contribute to p in a single raw
  if(other){
    results_b=results[1:nsp,]
    rBVI=colSums(results[(nsp+1):dim(results)[1], 2:dim(results)[2]])
    results_b[nsp+1,2:dim(results)[2]]=rBVI
    results_b$spp[nsp+1]="Others"
    results=results_b
    return(results)
  }
  
  #Return the results
  return(results)
}


