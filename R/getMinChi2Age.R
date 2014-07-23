findMinChi2AgeBB <- function(FTdataset){
  
  ages=FTdataset$grainAges
  nIList=FTdataset$nI
  nSList=FTdataset$nS
  
  grainAgeOrder = order(ages)
  nI = nIList[grainAgeOrder]
  nS = nSList[grainAgeOrder]

  N = length(nI)
  
  df = seq(0,N-1)
  Chi2_List=c()
  pValue_List=c()
  
  for(nGrains in 1:N){
  
  nImean = mean(nI[1:nGrains])
  nSmean = mean(nS[1:nGrains])
  
  Chi2=0
  for(i in 1:length(nGrains)){
    Chi2 = Chi2 + (nS[i]*nImean - nI[i]*nSmean)**2 / (nI[i]+nS[i]) 
  }
  
  Chi2_List[nGrains] = Chi2/(nImean*nSmean)
  pValue_List[nGrains] = pchisq(Chi2, df[nGrains], lower.tail=FALSE)
}
  
  
#   if(pValue<0.01){
#     comment="pValue<0.01 : Strong evidence against a common age model"
#   } else if (pValue<0.05) {
#     comment="pValue<0.05 : Moderate evidence against a common age model"
#   } else {
#     comment = "pValue>0.05 : Data consistent with a common age model"
#   }
#   
#   tmp = list(Chi2=Chi2 , df=df , pValue=pValue, significant=comment)
#   return(tmp)
  

  
  result=data.frame(ages[grainAgeOrder],nI,nS,Chi2_List, pValue_List, df)
  
  return(result)
}


findMinChi2Age <- function(FTdataset){
  
  ages=FTdataset$grainAges
  nIList=FTdataset$nI
  nSList=FTdataset$nS
  
  grainAgeOrder = order(ages)
  nI_All = nIList[grainAgeOrder]
  nS_All = nSList[grainAgeOrder]
  
  N = length(nI)
  
  df = seq(0,N-1)
  Chi2_List=c()
  pValue_List=c()
  pooledAge=c()
  
  zetaStErr = FTdataset$Zeta * FTdataset$relErrZeta
  
  for(nGrains in 1:N){
    nI = nI_All[1:nGrains]
    nS = nS_All[1:nGrains]
    
    nImean = mean(nI)
    nSmean = mean(nS)
    
    Chi2=0
    for(i in 1:length(nI)){
      Chi2 = Chi2 + (nS[i]*nImean - nI[i]*nSmean)**2 / (nI[i]+nS[i]) 
    }
    Chi2 = Chi2/(nImean*nSmean)
    pValue = pchisq(Chi2, df[nGrains], lower.tail=FALSE)
    
    Chi2_List[nGrains] = Chi2
    pValue_List[nGrains] = pValue
    
    tmp = calcPooledAge(nS[1:nGrains], nI[1:nGrains], FTdataset$Zeta, zetaStErr, FTdataset$rhoD, FTdataset$nD)
    pooledAge[nGrains] = tmp$pooledAge
  }
  
  #   if(pValue<0.01){
  #     comment="pValue<0.01 : Strong evidence against a common age model"
  #   } else if (pValue<0.05) {
  #     comment="pValue<0.05 : Moderate evidence against a common age model"
  #   } else {
  #     comment = "pValue>0.05 : Data consistent with a common age model"
  #   }
  #   
  #   tmp = list(Chi2=Chi2 , df=df , pValue=pValue, significant=comment)
  #   return(tmp)
   isCommonMinimum = (pValue_List>0.05)
  
  
  result=data.frame(ages[grainAgeOrder],nI_All,nS_All,Chi2_List, pValue_List, df,pooledAge,isCommonMinimum)
  
  return(result)
}