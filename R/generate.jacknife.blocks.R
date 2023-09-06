#' Generate block coordinates for block-jackknife
#' @param x A pooldata or countdata object containing SNP positions (snp.info slot)
#' @param nsnp.per.bjack.block Number of consecutive SNPs of each block-jackknife block 
#' @param verbose If TRUE extra information is printed on the terminal
#' @return A list with the two following elements:
#' \enumerate{
#' \item "blocks.det": A matrix with three columns containing for each identified block (in row) the index of the start SNP, the index of the end SNP and the block Size in bp
#' \item "snp.block.id": A vector containing the blocks assigned to each SNP eligible for block-Jacknife (non eligible SNPs ares assigned NA)
#' \item "nblocks": A scalar corresponding to the number of blocks
#' \item "nsnps": Number of SNPs eligible for block-jackknife 'i.e., included in one block
#' }
#' @examples
#' make.example.files(writing.dir=tempdir())
#' pooldata=popsync2pooldata(sync.file=paste0(tempdir(),"/ex.sync.gz"),poolsizes=rep(50,15))
#' bjack.blocks=generate.jackknife.blocks(pooldata,nsnp.per.bjack.block=50)
#' @export
generate.jackknife.blocks<-function(x,nsnp.per.bjack.block,verbose=TRUE){
  tmp.pos=x@snp.info[,2]
  det.idx.per.chr=matrix(unlist(by(1:x@nsnp,x@snp.info[,1],range)),ncol=2,byrow=T)
  if(nrow(det.idx.per.chr)==0){#in case all contig names are NA
    stop("Exit function: No chr/contigs available (information on SNP contig name might not have been provided)")
  }
  tmp.diff.pos=by(tmp.pos,x@snp.info[,1],diff,na.rm=T)
  det.idx.per.chr=cbind(det.idx.per.chr,
                        det.idx.per.chr[,2]-det.idx.per.chr[,1]+1,
                        tmp.pos[det.idx.per.chr[,2]]-tmp.pos[det.idx.per.chr[,1]],
                        matrix(unlist(lapply(tmp.diff.pos,f<-function(x){if(length(x)>0){return(range(x))}else{return(rep(NA,2))}})),ncol=2,byrow=T))
  colnames(det.idx.per.chr)=c("snp_idx1","snp_idx2","nsnps","size","min_intersnp_dist","max_intersnp_dist")
  tmp=sum(det.idx.per.chr[,5]<0,na.rm=T)
  if(tmp>0){
    warning(paste(tmp,"contigs with unordered SNPs positions (these will be discarded)\n"))
  }
  ctg.sel=det.idx.per.chr[,3]>nsnp.per.bjack.block & det.idx.per.chr[,5]>=0 #last=remove contigs unordered
  if(sum(ctg.sel)==0){
    stop("Exit function: No contig available after applying filtering steps (e.g., try lowering nsnp.per.bjack.block)\n")
  }
  det.idx.per.chr=det.idx.per.chr[ctg.sel,]
  if(sum(ctg.sel)==1){
#    ref.snp.idx=det.idx.per.chr[,1]:(x[2]-nsnp.per.bjack.block+1)
    block.start=seq(det.idx.per.chr[1],det.idx.per.chr[2]-nsnp.per.bjack.block+1,nsnp.per.bjack.block)
#    snp.idx.to.sample=det.idx.per.chr[ctg.sel,1]:(det.idx.per.chr[ctg.sel,2]-nsnp.per.bjack.block+1)
  }else{
#    ref.snp.idx=unlist(apply(det.idx.per.chr[,1:2],1,f<-function(x){return(x[1]:(x[2]-nsnp.per.bjack.block+1))})) #we only keep to compute stats those snps that could be include in a block (not to bias estimates if a lot of small ctgs)
    block.start=as.vector(unlist(apply(det.idx.per.chr[,1:2],1,f<-function(x){return(seq(x[1],x[2]-nsnp.per.bjack.block+1,nsnp.per.bjack.block))}))) 
#    snp.idx.to.sample=unlist(apply(det.idx.per.chr[,1:2],1,f<-function(x){return(x[1]:(x[2]-nsnp.per.bjack.block+1))}))
  }
  nblocks=length(block.start)
  if(nblocks<=1){
    stop("Exit function: <=1 block available\n")
  }
  blocks.det=cbind(block.start,block.start+nsnp.per.bjack.block-1)  
  blocks.det=cbind(blocks.det,tmp.pos[blocks.det[,2]]-tmp.pos[blocks.det[,1]])
  colnames(blocks.det)=c("idx.start","idx.end","size")
  ##vector of block id
  snp.block.id=rep(NA,x@nsnp)
  for(i in 1:nblocks){snp.block.id[blocks.det[i,1]:blocks.det[i,2]]=i}
  if(verbose){
   cat(nblocks,"Jackknife blocks identified with",sum(!is.na(snp.block.id)),"SNPs (out of",x@nsnp,").\n SNPs map to",nrow(det.idx.per.chr),"different chrom/scaffolds\n")
   cat("Average (min-max) Block Sizes:",round(mean(blocks.det[,3]*1e-6),3),"(",round(min(blocks.det[,3]*1e-6),3),"-",round(max(blocks.det[,3]*1e-6),3),") Mb\n")
  }
  return(list(blocks.det=blocks.det,snp.block.id=snp.block.id,nblocks=nblocks,nsnps=sum(!is.na(snp.block.id))))
}