This is the guidance of PRS construction using transfer learning.
It can help predict disease risk in minor ancestry better using summary statistics from other resources.

The main function "PRS_TransferLearning" contains the following inputs:
 
ped_file: This file contains the information of FID, IID, Y and covariates.
Covar_name:The name of covariates we need to adjust in the model, such as c("Sex","BY","P1","P2","P3","P4"). Note that all names must corrspond to the columns in the ped file.
Y_name: The name of Y in the model, such as "LDL". Note that the Y name must corrspond to a column in the ped file.
Ytype: The type of Y should be either "C"(continuous) or "B"(binary).
train_file:The prefix of plink file of the training data in the target population. Note that we use the training data to train the new effect sizes of the target population. 
test_file: The prefix of plink file of the testing data in the target population. Note that we use the testing data to choose the best tuning parameter. The user should use the other independent data for validation. 
sum_stats_file:The file contains the effect sizes from other resources. Specifically it contains the following three columns:"SNP","A1","Beta". "SNP" is the SNPID (the same format as SNPID in plink files); "A1" is the alternative (effect) allele; "Beta" is the effect size. 
Note that the effect size in the file is different from effect sizes from GWAS. The effect size here should come from the output of PRS methods, such as lassosum and PRSCS.
LDblocks: This will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome. Currently we support three types:"EUR.hg19","AFR.hg19","ASN.hg19", corresponding to European, African and Asian populations, respectively.
tempfile: The prefix of the file location which can be used to store the temporary files. Note that the function needs to save files in this directory.

The current output contains four components:

best.learning.rate: the learning rate we can use in order to achieve the best risk prediction.
best.iteration: the number of iterations we should stop in order to achieve the best risk prediction.
best.beta: the data frame containing three columns: "SNP","A1","beta". Note that this is the best effect size we can use to construct PRS, selected using best.learning.rate and best.iteration.  
best.PRS: This component provides PRS for testing file. It is a data frame containing four columns:"FID","IID","PRS.NULL","PRS.TL". Note that "PRS.NULL" is calculated based on effect sizes provided by sum_stats_file and "PRS.TL" is calculated based on best.beta. 
param_table: the data frame containing a grid of candidates of learning rates and the number of iterations that we consider. 


Here is one test example:
##plink-1.9 needs to be pre-installed.
source("/net/snowwhite/home/zczhao/PRS/Pipeline/Rpackage/PRS_TransferLearning.R")
library(data.table)
library(lassosum)
library(Matrix)
library(parallel)
ped_file="/net/snowwhite/home/zczhao/PRS/Pipeline/Pheno/AFR_8traits_1113.ped";
Covar_name=c("Sex","BY","P1","P2","P3","P4");Y_name="LDL";Ytype="C"
train_file="/net/snowwhite/home/zczhao/PRS/Pipeline/data/African_4kGWAS_plink"
test_file="/net/csgspare3/snowwhite.archive/zczhao/PRS_Geno/African/African_2ktrain_plink"
sum_stats_file=paste0("/net/csgspare3/snowwhite.archive/zczhao/PRS_Geno/Method1/Training/out_","lassosum","_","AFR_5k","_","LDL","_","eur",".txt")
LDblocks="EUR.hg19"
tempfile="/net/snowwhite/home/zczhao/PRS/Pipeline/Rpackage/Temp/LDL_0408_"
system.time({out.beta=PRS_TransferLearning(ped_file,Covar_name,Y_name,Ytype, train_file,test_file,sum_stats_file,LDblocks,tempfile)})
summary(out.beta)



