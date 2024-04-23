# TL-PRS
This R package helps users to construct multi-ethnic polygenic risk score (PRS) using transfer learning. It can help predict PRS of minor ancestry using summary statistics from exsiting resources, such as UK Biobank.

This package contains two main function: TL_PRS and TL_PRS_Ind. The main difference between the two functions is that TL_PRS_Ind requires individual-level data for training while TL_PRS doen't. But both two functions require individual-level data for validation. According to our simulation results, TL_PRS is more recommended becasue it doesn't requirement individual-level data for training and it has similar performance as TL_PRS_Ind. Here we show how to implment TL_PRS in R software.

## Reference
[Zhao, Z., Fritsche, L.G., Smith, J.A., Mukherjee, B. and Lee, S., 2022. The Construction of Multi-ethnic Polygenic Risk Score using Transfer Learning. medRxiv.](https://www.medrxiv.org/content/10.1101/2022.03.08.22272114v1)

## Installation
`TL-PRS` requires the software 'plink' as well as the following R packages:  `lassosum` and `parallel`. Install them by: 

```r
install.packages(c("lassosum", "parallel"), dependencies=TRUE)
```

If you have `devtools`, you can type: 
```r
install_github("ZhangchenZhao/TLPRS")
```
for the latest development version. Or you can clone the latest development version here and install yourself using `devtools`. 

## Inputs of TL_PRS
1. `ped_file`:
The location path of ped file, where contains the information of FID, IID, outcome (Y) and covariates. Note that the file only requires samples from `validate_file`.

2. `Covar_name`:
A vector of names of covariates we need to adjust in the model, such as c("Sex","Age"). Note that all names must corrspond to the columns in the ped file.

3. `Y_name`: 
The name of Y in the model, such as "LDL". Note that the Y name must corrspond to a column in the ped file.

4. `Ytype`: 
The type of Y should be either "C"(continuous) or "B"(binary).

5. `train_file`:
The prefix of plink file of the training data in the target population. Note that we use the training data to train the new effect sizes of the target population. 

6. `validate_file`:
The prefix of plink file of the validation data in the target population. Note that we use the validation data to choose the best tuning parameter. 

7. `sum_stats_file`:
The location path of effect size file. We usually can obtain this file by existing PRS methods, such as lassosum/PRS-CS. Specifically it contains the following three columns:"SNP","A1","Beta". "SNP" is the SNPID (the same format as SNPID in plink files); "A1" is the alternative (effect) allele; "Beta" is the effect size. 

8. `target_sumstats_file`:
The location path of summary stats file for the target population. The file requires the following columns: "SNP", "A1", "beta", "N", and "p", wheren "beta" is the effect size from summary statistics, "N" is the sample size used for calculating summary statistics, and "p" is p-value of the SNP. 

9. `LDblocks`:
This will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome. Currently we support three types:"EUR.hg19","AFR.hg19","ASN.hg19", corresponding to European, African and Asian populations, respectively.

10. `outfile`:
The prefix of the file location which can be used to store the final output files. Note that the function needs to save files in this directory.

11. `cluster`:
A cluster object from the parallel package for parallel computing

## Outputs of TL_PRS
1. `best.learning.rate`: 
the learning rate we can use in order to achieve the best risk prediction.

2. `best.iteration`: 
the number of iterations we should stop in order to achieve the best risk prediction.

3. `best.beta`: 
the data frame containing three columns: "SNP","A1","beta". Note that this is the best effect size we can use to construct PRS, selected using best.learning.rate and best.iteration.  

4. `best.PRS`: 
This component provides PRS for testing file. It is a data frame containing four columns:"FID","IID","PRS.NULL","PRS.TL". Note that "PRS.NULL" is calculated based on effect sizes provided by sum_stats_file and "PRS.TL" is calculated based on best.beta. 

5. `param_table`: 
the data frame containing a grid of candidates of learning rates and the number of iterations that we consider. 


## Example of TL_PRS 

I list one example script below. If you would like to play with a toy example and run TL-PRS using this example script, please use the link to download the data:
https://www.dropbox.com/sh/40vewd1kuxcbeev/AAD7Dj3H-sBTWv2ObUIDEHFya?dl=0 

```
##plink-1.9 needs to be pre-installed.

#########Step 1. Run TL-PRS using example data#########
library(data.table)
library(lassosum)
library(TLPRS)
setwd("XXX/XXX")  ###setup your work path.
ped_file="ped_validation.txt";
Covar_name="sex"
Y_name="Y"
Ytype="C"
train_file="ExampleData_1000G_African_train"
validate_file="ExampleData_1000G_African_validate"
sum_stats_file="Beta_lassosum_usingEuropeanSummaryStatistics.txt"
target_sumstats_file="African_SummaryStatistics.txt"
LDblocks="AFR.hg19"
outfile="Output_TLPRS"
system.time({out.beta=TL_PRS(ped_file,Covar_name,Y_name, Ytype="C",train_file,validate_file,sum_stats_file,target_sumstats_file, LDblocks,outfile,cluster=NULL)})
summary(out.beta)
write.table(out.beta, file=,quote=F,row.names=F,col.names=F)

#########Step 2. Compare TL-PRS with its baseline method using example testing data##########
test_file="ExampleData_1000G_African_test"
ped_test_file="ped_test.txt"

Correlation_betweenPRSandY <- function( test_file, stats_file, ped_test_file,Y_name){
	Y_test=read.table(ped_test_file,header=T)
	cmd=paste0("plink-1.9 --bfile ",test_file,"  --score ",stats_file, " sum", " --out ",stats_file,".test.PRS")
	system(cmd)
	temp=read.table(paste0(stats_file,".test.PRS.profile"),header=T)
	merged=merge(temp,Y_test,by=c("FID","IID"))
	return(cor(merged[,Y_name],merged$SCORESUM))
}

###Correlation between lassosum PRS and Y in the testing data
Correlation_betweenPRSandY( test_file,sum_stats_file, ped_test_file, Y_name)  

###Correlation between TL-PRS and Y in the testing data
Correlation_betweenPRSandY( test_file, paste0(outfile, "_best.beta.txt"), ped_test_file,Y_name)

```

## Support
If there are any further questions or problems with running or installing `TL-PRS`, please do email me at <zczhao@umich.edu>. 
