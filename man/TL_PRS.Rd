\name{TL_PRS}
\alias{TL_PRS}
\title{
Perform TL-PRS using summary statistics
}
\description{
  Perform TL-PRS using summary statistics
}
\usage{
TL_PRS(ped_file,Covar_name,Y_name, Ytype="C",train_file,test_file,sum_stats_file,
       target_sumstats_file, LDblocks="EUR.hg19",outfile,cluster=NULL)
}
\arguments{
  \item{ped_file}{
The location path of ped file, where contains the information of FID, IID, outcome (Y) and covariates. Note that the file only requires samples from test_file.
}
  \item{Covar_name}{
A vector of names of covariates we need to adjust in the model, such as c("Sex","BY","PC1","PC2","PC3","PC4"). Note that all names must corrspond to the columns in the ped file.
}
  \item{Y_name}{
The name of Y in the model, such as "LDL". Note that the Y name must corrspond to a column in the ped file.
}
  \item{Ytype}{
The type of Y should be either "C"(continuous) or "B"(binary).
}
  \item{train_file}{
The prefix of plink file of the training data in the target population. Note that we use the training data to train the new effect sizes of the target population. 
}
  \item{test_file}{
The prefix of plink file of the testing data in the target population. Note that we use the testing data to choose the best tuning parameter. 
}
  \item{sum_stats_file}{
The location path of effect size file. We usually can obtain this file by existing PRS methods, such as lassosum/PRS-CS. Specifically it contains the following three columns:"SNP","A1","Beta". "SNP" is the SNPID (the same format as SNPID in plink files); "A1" is the alternative (effect) allele; "Beta" is the effect size. 
}
  \item{target_sumstats_file}{
The location path of summary stats file for the target population. The file requires the following columns: "SNP", "A1", "beta", "N", and "p", wheren "beta" is the effect size from summary statistics, "N" is the sample size used for calculating summary statistics, and "p" is p-value of the SNP. 
}
  \item{LDblocks}{
This will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome. Currently we support three types:"EUR.hg19","AFR.hg19","ASN.hg19", corresponding to European, African and Asian populations, respectively.
}
  \item{outfile}{
The prefix of the file location which can be used to store the final output files. Note that the function needs to save files in this directory.
}
  \item{cluster}{
A \code{cluster} object from the \code{parallel} package for parallel computing
}
}

\details{
\code{genos} can have discrete \code{0, 1, 2} values or continuous values between \code{[0,2]}. The genotype or dosage values can represent any of the major allele, minor allele, reference allele or alternate allele counts (or dosages), as long as it is consistent throughout the subjects.

\code{genos} can have missing values denoted by the \code{missing.id} argument. Such missing values will be imputed using mean imputation. \code{pheno} or \code{cov} cannot have missing values. 

\code{pheno} and \code{cov} are ignored if \code{obj.null} is provided. If both \code{obj.null} and \code{cov} is missing, or \code{obj.null} is missing and \code{cov=NULL}, then the vector \code{rep(1,n)} is assigned to \code{cov}, where \code{n} is the number of subjects.

\code{method} = "SPA" is the basic saddlepoint approximation based test without the partially normal approximation improvement.
\code{method} = "fastSPA" utilizes the partially normal approximation approach for improved efficiency, especially for rare variants.

Beta parameters are estimated using Firth's method, and thus computationally expensive. Therefore, it is recommended that beta parameters are only to be estimated when the p-value is very small (denoted by \code{beta.Cutoff}). The code for beta estimation is as implemented by Clement Ma in the EPACTS software.
}

\value{
\item{best.learning.rate}{the learning rate we can use in order to achieve the best risk prediction.}
\item{best.iteration}{the number of iterations we should stop in order to achieve the best risk prediction.}
\item{best.beta}{the data frame containing three columns: "SNP","A1","beta". Note that this is the best effect size we can use to construct PRS, selected using best.learning.rate and best.iteration.  
}
\item{best.PRS}{This component provides PRS for testing file. It is a data frame containing four columns:"FID","IID","PRS.NULL","PRS.TL". Note that "PRS.NULL" is calculated based on effect sizes provided by sum_stats_file and "PRS.TL" is calculated based on best.beta. 
}
\item{param_table}{the data frame containing a grid of candidates of learning rates and the number of iterations that we consider. }
}
\author{
Zhangchen Zhao, \email{zczhao@umich.edu}
}
\references{
Zhao, Z., Fritsche, L.G., Smith, J.A., Mukherjee, B. and Lee, S., 2022. The Construction of Multi-ethnic Polygenic Risk Score using Transfer Learning. medRxiv.
}


