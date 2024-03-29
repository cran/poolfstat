% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf2pooldata.R
\name{vcf2pooldata}
\alias{vcf2pooldata}
\title{Convert a VCF file into a pooldata object.}
\usage{
vcf2pooldata(
  vcf.file = "",
  poolsizes = NA,
  poolnames = NA,
  min.cov.per.pool = -1,
  min.rc = 1,
  max.cov.per.pool = 1e+06,
  min.maf = -1,
  remove.indels = FALSE,
  min.dist.from.indels = 0,
  nlines.per.readblock = 1e+06,
  verbose = TRUE
)
}
\arguments{
\item{vcf.file}{The name (or a path) of the Popoolation sync file (might be in compressed format)}

\item{poolsizes}{A numeric vector with haploid pool sizes}

\item{poolnames}{A character vector with the names of pool}

\item{min.cov.per.pool}{Minimal allowed read count (per pool). If at least one pool is not covered by at least min.cov.perpool reads, the position is discarded}

\item{min.rc}{Minimal allowed read count per base (options silenced for VarScan vcf). Bases covered by less than min.rc reads are discarded and considered as sequencing error. For instance, if nucleotides A, C, G and T are covered by respectively 100, 15, 0 and 1 over all the pools, setting min.rc to 0 will lead to discard the position (the polymorphism being considered as tri-allelic), while setting min.rc to 1 (or 2, 3..14) will make the position be considered as a SNP with two alleles A and C (the only read for allele T being disregarded). For VarScan vcf, markers with more than one alternative allele are discarded because the VarScan AD field only contains one alternate read count.}

\item{max.cov.per.pool}{Maximal allowed read count (per pool). If at least one pool is covered by more than min.cov.perpool reads, the position is discarded}

\item{min.maf}{Minimal allowed Minor Allele Frequency (computed from the ratio overall read counts for the reference allele over the read coverage)}

\item{remove.indels}{Remove indels identified using the number of characters of the alleles in the REF or ALT fields (i.e., if at least one allele is more than 1 character, the position is discarded)}

\item{min.dist.from.indels}{Remove SNPs within min.dist.from.indels from an indel i.e. SNP with position p verifying (indel.pos-min.dist)<=p<=(indel.pos+min.dist+l.indels-1) where l.indel=length of the ref. indel allele. If min.dist.from.indels>0, INDELS are also removed (i.e., remove.indels is set to TRUE).}

\item{nlines.per.readblock}{Number of Lines read simultaneously. Should be adapted to the available RAM.}

\item{verbose}{If TRUE extra information is printed on the terminal}
}
\value{
A pooldata object containing 7 elements:
\enumerate{
\item "refallele.readcount": a matrix with nsnp rows and npools columns containing read counts for the reference allele (chosen arbitrarily) in each pool
\item "readcoverage": a matrix with nsnp rows and npools columns containing read coverage in each pool
\item "snp.info": a matrix with nsnp rows and four columns containing respectively the contig (or chromosome) name (1st column) and position (2nd column) of the SNP; the allele taken as reference in the refallele.readcount matrix (3rd column); and the alternative allele (4th column)
\item "poolsizes": a vector of length npools containing the haploid pool sizes
\item "poolnames": a vector of length npools containing the names of the pools
\item "nsnp": a scalar corresponding to the number of SNPs
\item "npools": a scalar corresponding to the number of pools
}
}
\description{
Convert VCF files into a pooldata object.
}
\details{
Genotype format in the vcf file for each pool is assumed to contain either i) an AD field containing allele counts separated by a comma (as produced by popular software such as GATK or samtools/bcftools) or ii) both a RD (reference allele count) and a AD (alternate allele count) as obtained with the VarScan mpileup2snp program (when run with the --output-vcf option). The underlying format is automatically detected by the function. For VarScan generated vcf, it should be noticed that SNPs with more than one alternate allele are discarded (because only a single count is then reported in the AD fields) making the min.rc unavailable. The VarScan --min-reads2 option might replace to some extent this functionality although SNP where the two major alleles in the Pool-Seq data are different from the reference allele (e.g., expected to be more frequent when using a distantly related reference genome for mapping) will be disregarded.
}
\examples{
 make.example.files(writing.dir=tempdir())
 pooldata=vcf2pooldata(vcf.file=paste0(tempdir(),"/ex.vcf.gz"),poolsizes=rep(50,15))
}
