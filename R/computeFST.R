#' Compute FST from Pool-Seq data
#' @param pooldata A pooldata object containing Pool-Seq information
#' @param method Either "Anova" (default method as described in the manuscript) or "Identity" (relies on an alternative modeling consisting in estimating unbiased Probability of Identity within and across pairs of pools)
#' @param snp.index A list of SNP to be considered in the computation (by default all the SNP are considered)
#' @return A list with the four following elements:
#' \enumerate{
#' \item "FST": a scalar corresponding to the estimate the global FST
#' \item "snp.FST": a vector containing estimates of SNP-specific FST
#' \item "snp.Q1": a vector containing estimates of the overall within pop. SNP-specific probability of identity
#' \item "snp.Q2": a vector containing estimates of the overall between pop. SNP-specific probability of identity
#' }
#' @seealso To generate pooldata object, see \code{\link{vcf2pooldata}}, \code{\link{popsync2pooldata}}
#' @examples
#'  make.example.files(writing.dir=tempdir())
#'  pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#'  res.fst=computeFST(pooldata)
#' @export
computeFST=function(pooldata,method="Anova",snp.index=NA){
  if(!(method %in% c("Identity","Anova"))){stop("method should either be Identity or Anova (default)")}
	if(!(is.pooldata(pooldata))) {stop("The data are not formatted as a valid pooldata object...
	  (see the popsync2pooldata, vcf2pooldata, genobaypass2pooldata and selestim2pooldata functions)")} 
  npop=pooldata@npools ; nsnp=pooldata@nsnp
  if(is.na(snp.index)){snp.index=1:nsnp}
  YY=pooldata@refallele.readcount[snp.index,]
  NN=pooldata@readcoverage[snp.index,]
  poolsize=pooldata@poolsizes

 if(method=="Identity"){
  Q1=as.matrix((YY*(YY-1) + (NN-YY)*(NN-YY-1) )/(NN*(NN-1)))
  Q1 = (1/(matrix(1,nsnp,npop) %*% diag(poolsize-1)))*(Q1 %*% diag(poolsize) - 1)
  lambdaj=poolsize*(poolsize-1)
  lambdaj=lambdaj/sum(lambdaj)
  snp.Q1=rowSums(Q1%*%diag(lambdaj))
  hat.Q1=mean(snp.Q1,na.rm=T) 

  Q2=matrix(0,nrow(YY),npop*(npop-1)/2)
  omegajj=rep(0,npop*(npop-1)/2)
  cnt=0
  for(i in 1:(npop-1)){
    for(j in (i+1):npop){
      cnt=cnt+1
      omegajj[cnt]=poolsize[i]*poolsize[j]
      Q2[,cnt]=(YY[,i]*YY[,j] + (NN[,i]-YY[,i])*(NN[,j]-YY[,j]))/(NN[,i]*NN[,j])
    }
  }
  snp.Q2=rowSums(Q2%*%diag(omegajj/sum(omegajj)))
  hat.Q2=mean(snp.Q2,na.rm=TRUE)
  
  rslt=list(FST=(hat.Q1-hat.Q2)/(1-hat.Q2),snp.Q1=snp.Q1,snp.Q2=snp.Q2,snp.FST=(snp.Q1-snp.Q2)/(1-snp.Q2) )
 }

	if (method=="Anova"){
		mtrx.n_i <- matrix(poolsize,nrow = nsnp,ncol = npop,byrow = TRUE)
		C_1 <- rowSums(NN)
		C_2 <- rowSums(NN^2)
		D_2 <- rowSums(NN / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i,na.rm = TRUE)
		D_2.star <- rowSums(NN * (NN / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i),na.rm = TRUE) / C_1
		n_c <- (C_1 - C_2 / C_1) / (D_2 - D_2.star)
		YY_alt <- NN - YY
		SSI <- rowSums(YY - YY^2 / NN + YY_alt - YY_alt^2 / NN,na.rm = TRUE)
		SSP <- rowSums(NN * ((YY / NN) - (rowSums(YY) / C_1))^2 + NN * ((YY_alt / NN) - (rowSums(YY_alt) / C_1))^2,na.rm = TRUE)
		MSI <- SSI / (C_1 - D_2)
		MSP <- SSP / (D_2 - D_2.star)
		F_ST <- (MSP - MSI)  / (MSP + (n_c - 1) * MSI)
		keep <- !is.na(F_ST)
		F_ST_multi <- sum(MSP[keep] - MSI[keep])  / sum(MSP[keep] + (n_c[keep] - 1) * MSI[keep])
		Q_1 <- 1 - MSI
		Q_2 <- 1 - MSI - (MSP - MSI) / n_c
		rslt <- list(snp.FST = F_ST,snp.Q1 = Q_1,snp.Q2 = Q_2,FST = F_ST_multi)
	}
	return(rslt)
}
