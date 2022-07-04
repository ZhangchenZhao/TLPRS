library(data.table)
library(lassosum)
library(Matrix)
library(parallel)


##calculate nagelkerke R2 for binary phenotypes
nagelkerke<-function (fit, null = NULL, restrictNobs = FALSE){
    TOGGLE = (class(fit)[1] == "lm" | class(fit)[1] == "gls" |
        class(fit)[1] == "lme" | class(fit)[1] == "glm" | class(fit)[1] ==
        "negbin" | class(fit)[1] == "zeroinfl" | class(fit)[1] ==
        "clm" | class(fit)[1] == "vglm" | class(fit)[1] == "betareg" |
        class(fit)[1] == "rq")
    BOGGLE = (class(fit)[1] == "nls" | class(fit)[1] == "lmerMod" |
        class(fit)[1] == "glmerMod" | class(fit)[1] == "merModLmerTest" |
        class(fit)[1] == "lmerModLmerTest" | class(fit)[1] ==
        "clmm")
    SMOGGLE = (class(fit)[1] == "lmerMod" | class(fit)[1] ==
        "glmerMod" | class(fit)[1] == "merModLmerTest" | class(fit)[1] ==
        "lmerModLmerTest" | class(fit)[1] == "vglm")
    ZOGGLE = (class(fit)[1] == "zeroinfl")
    ZOGGLE2 = (class(fit)[1] == "rq")
    NOGGLE = is.null(null)
    ERROR = "Note: For models fit with REML, these statistics are based on refitting with ML"
    ERROR2 = "None"
    if (!restrictNobs & NOGGLE & TOGGLE) {
        null = update(fit, ~1)
    }
    if (restrictNobs & NOGGLE & TOGGLE) {
        null = update(fit, ~1, data = fit$model)
    }
    if (restrictNobs & !NOGGLE) {
        null = update(null, data = fit$model)
    }
    if (NOGGLE & BOGGLE) {
        ERROR = "You need to supply a null model for nls, lmer, glmer, or clmm"
    }
    if ((!TOGGLE) & (!BOGGLE)) {
        ERROR = "This function will work with lm, gls, lme, lmer, glmer, glm, negbin, zeroinfl, nls, clm, clmm, and vglm"
    }
    SMOGGLE2 = (class(null)[1] == "lmerMod" | class(null)[1] ==
        "glmerMod" | class(null)[1] == "merModLmerTest" | class(null)[1] ==
        "lmerModLmerTest" | class(null)[1] == "vglm")
    Y = matrix(rep(NA, 2), ncol = 1)
    colnames(Y) = ""
    rownames(Y) = c("Model:", "Null:")
    Z = matrix(rep(NA, 3), ncol = 1)
    colnames(Z) = c("Pseudo.R.squared")
    rownames(Z) = c("McFadden", "Cox and Snell (ML)", "Nagelkerke (Cragg and Uhler)")
    X = matrix(rep(NA, 4), ncol = 4)
    colnames(X) = c("Df.diff", "LogLik.diff", "Chisq", "p.value")
    rownames(X) = ""
    U = matrix(rep(NA, 2), ncol = 1)
    colnames(U) = ""
    rownames(U) = c("Model:", "Null:")
    if (TOGGLE | BOGGLE) {
        if (!SMOGGLE) {
            Y[1] = toString(fit$call)
        }
        if (SMOGGLE) {
            Y[1] = toString(fit@call)
        }
    }
    if (TOGGLE | (BOGGLE & !NOGGLE)) {
        if (!SMOGGLE2) {
            Y[2] = toString(null$call)
        }
        if (SMOGGLE2) {
            Y[2] = toString(null@call)
        }
        if (!ZOGGLE & !ZOGGLE2) {
            N = nobs(fit)
            U[1, 1] = nobs(fit)
            U[2, 1] = nobs(null)
        }
        if (!ZOGGLE & ZOGGLE2) {
            N = length(fit$y)
            U[1, 1] = length(fit$y)
            U[2, 1] = length(null$y)
        }
        if (ZOGGLE) {
            N = fit$n
            U[1, 1] = fit$n
            U[2, 1] = null$n
        }
        if (U[1, 1] != U[2, 1]) {
            ERROR2 = "WARNING: Fitted and null models have different numbers of observations"
        }
        m = suppressWarnings(logLik(fit, REML = FALSE))[1]
        n = suppressWarnings(logLik(null, REML = FALSE))[1]
        mf = 1 - m/n
        Z[1, ] = signif(mf, digits = 6)
        cs = 1 - exp(-2/N * (m - n))
        Z[2, ] = signif(cs, digits = 6)
        nk = cs/(1 - exp(2/N * n))
        Z[3, ] = signif(nk, digits = 6)
        o = n - m
        dfm = attr(logLik(fit), "df")
        dfn = attr(logLik(null), "df")
        if (class(fit)[1] == "vglm") {
            dfm = df.residual(fit)
        }
        if (class(fit)[1] == "vglm") {
            dfn = df.residual(null)
        }
        dff = dfn - dfm
        CHI = 2 * (m - n)
        P = pchisq(CHI, abs(dff), lower.tail = FALSE)
        X[1, 1] = dff
        X[1, 2] = signif(o, digits = 5)
        X[1, 3] = signif(CHI, digits = 5)
        X[1, 4] = signif(P, digits = 5)
    }
    W = ERROR
    WW = ERROR2
    V = list(Y, Z, X, U, W, WW)
    names(V) = c("Models", "Pseudo.R.squared.for.model.vs.null",
        "Likelihood.ratio.test", "Number.of.observations", "Messages",
        "Warnings")
    return(V)
}




##resulttemp,ped,Covar_name,Y_name
logistic_result_generator_old<-function(PRS,ped,Covar_name,Y_name){
        datatemp2=merge(PRS,ped,by.x=c("FID","IID"),by.y=c("FID","IID"))
  	datatemp2$PRS=scale(datatemp2$SCORESUM, center = TRUE, scale = TRUE)
	args=gsub(",","+",toString(Covar_name))
        modeltemp0=glm(paste0(Y_name,"~",args), data=datatemp2, family = "binomial")
        modeltemp=glm(paste0(Y_name,"~ PRS +",args), data=datatemp2, family = "binomial")
	#casecontrol=paste0(sum(datatemp2$pheno==1,na.rm=T),":",sum(datatemp2$pheno==0,na.rm=T))
	r2=1-modeltemp$deviance/modeltemp0$deviance
	#ci=try(exp(cbind(coef(modeltemp), confint(modeltemp))["PRS",1:3] ),silent=TRUE)
	#if (class(ci)!="try-error"){
	#	or=ci[1]
	#	ci2=paste0("(",round(ci[2],2),",",round(ci[3],2),")")
	#} else {
	#	or=NA;ci2=NA
	#}
	#prs.p=coef(summary( modeltemp))["PRS",4]
	#
        #list.out=cbind(class=name1,casecontrol,r2,or,ci2, prs.p,r2)
        return(r2)
}

##Generate logitic regression results for binary phenotype
logistic_result_generator<-function(PRS,ped,Covar_name,Y_name){
        datatemp2=merge(PRS,ped,by.x=c("FID","IID"),by.y=c("FID","IID"))
  	datatemp2$PRS=scale(datatemp2$SCORESUM, center = TRUE, scale = TRUE)
	args=gsub(",","+",toString(Covar_name))
	if (args!=""){
        	modeltemp0=glm(paste0(Y_name,"~",args), data=datatemp2, family = "binomial")
        	modeltemp=glm(paste0(Y_name,"~ PRS +",args), data=datatemp2, family = "binomial")
	} else {
        	modeltemp0=glm(paste0(Y_name,"~1"), data=datatemp2, family = "binomial")
        	modeltemp=glm(paste0(Y_name,"~PRS"), data=datatemp2, family = "binomial")
	}
	r2=nagelkerke(modeltemp,null=modeltemp0)[[2]][3]
        return(r2)
}


##Generate linear regression results for continuous phenotype
linear_result_generator<-function(PRS,ped,Covar_name,Y_name){
	datatemp2=merge(PRS,ped,by.x=c("FID","IID"),by.y=c("FID","IID"))
	datatemp2$PRS=scale(datatemp2$SCORESUM, center = TRUE, scale = TRUE)
	args=gsub(",","+",toString(Covar_name))
	if (args!=""){
		modeltemp0=lm(paste0(Y_name,"~",args), data=datatemp2)
		modeltemp=lm(paste0(Y_name,"~ PRS +",args), data=datatemp2)
	} else{
		modeltemp0=lm(paste0(Y_name,"~1"), data=datatemp2)
		modeltemp=lm(paste0(Y_name,"~PRS"), data=datatemp2)
	}
	#prs.coef <- summary(modeltemp)$coeff[2,]
	#prs.beta <- as.numeric(prs.coef[1])
	#prs.se <- as.numeric(prs.coef[2])
	#prs.p <- as.numeric(prs.coef[4])
	model0.r2=summary(modeltemp0)$r.squared
	model1.r2=summary(modeltemp)$r.squared
	prs.r2 = summary(modeltemp)$r.squared-model0.r2
	#list.out=cbind(class=name1, prs.beta,prs.se,prs.p, model0.r2, model1.r2, prs.r2)
	return(prs.r2)
}

result_generator<-function(PRS,ped,Covar_name,Y_name,Ytype){
	if (Ytype=="C"){
		return(linear_result_generator(PRS,ped,Covar_name,Y_name))
	} else {
		return(logistic_result_generator(PRS,ped,Covar_name,Y_name))
	}
}

##test.bfile=train_file;beta.file=sum_stats_file;
##covar=Covar_name;kword=Y_name;temp.file=paste0(tempfile,"_step0") 
calculate_betaPRS<-function(test.bfile,beta.file,ped,covar,kword,temp.file){
	cmd=paste0("plink-1.9 --bfile ",test.bfile,"  --score ",beta.file, " sum", " --out ",temp.file,".train.PRS")
	system(cmd)
	temp=read.table(paste0(temp.file,".train.PRS.profile"),header=T)[,c(1,2,6)]
	temp2=merge(temp, ped,by.x=c("FID","IID"),by.y=c("FID","IID"),sort=F)
	list1=which(is.na(temp2[,kword]))
	if (length(list1)>0){
		temp3=temp2[-which(is.na(temp2[,kword])),]
	} else {temp3=temp2}
	temp3[,3]=temp3[,3]-mean(temp3[,3],na.rm=T)
	for (i in 4:ncol(temp3)){
		temp3[,i]=scale(temp3[,i],center=TRUE,scale=TRUE) ##Covariates and the phenotype
	}
	fit1=lm(paste0(kword,"~."),data=temp3[,-c(1,2)])

	covtemp=as.matrix(cbind(1,temp3[,-which(colnames(temp3)%in% c("FID","IID",kword,"SCORESUM"))]))
	
	fit1$RES_withoutPRS=temp3[,kword]-covtemp%*%fit1$coef[which(names(fit1$coef)!="SCORESUM")]


	return(list("model"=fit1,"data"=temp3))
}




Calculate_PRS<-function(test.bfile,B.beta.info,B.beta.all){
	test.bim=fread(paste0(test.bfile,".bim"),header=F)
	test.bim$V1 <- as.character(sub("^chr", "", test.bim$V1, ignore.case = T))
	test.bim$.index.ref <- 1:nrow(test.bim)

	B.beta.info$.index.tomatch <- 1:nrow(B.beta.info)

	merged <- merge(test.bim,   B.beta.info,all=F, 
                  by.x="V2", by.y="SNP",sort=F)
	###check the sign##
	list1=which(merged$A2==merged$V5)
	if (length(list1)>0){
		B.beta.all[list1,]=-B.beta.all[list1,]
	}

	flag<-rep(FALSE, nrow(test.bim))
	flag[merged$.index.ref]<-TRUE

	BM=as.matrix(B.beta.all[merged$.index.tomatch,])
	BM0=list()
	BM0[[1]]=BM

	system.time({
	pgs <- lapply(BM0, function(x) pgs(bfile=test.bfile, weights = x, 
           extract=flag, keep=NULL, 
           cluster=NULL))
	})

	fam=read.table(paste0(test.bfile,".fam"))[,c(1,2)]
	colnames(fam)[1:2]=c("FID","IID")
	out.all=cbind(fam,pgs[[1]])
	return(out.all)
}


PRStr_main_check<-function(ped_file,Covar_name,Y_name, Ytype,train_file,test_file,sum_stats_file,LDblocks){
	out1=0
	if (!file.exists(ped_file)){out1="The ped file does not exist!"} else {
		temp=fread(ped_file,header=T,nrow=1)
		if (!Y_name %in% colnames(temp)){out1="Y does not exist in the ped file!"}
		if (sum(!c("FID","IID") %in% colnames(temp))>0){out1="FID and IID do not exist in the ped file!"}
		if (sum(Covar_name=="")<1){
			if (sum(!Covar_name %in% colnames(temp))>0){out1="The covariates do not exist in the ped file!"}	
		}
	}
	if (!Ytype %in% c("C", "B")) {out1="The Y type is wrong! Now we only support continuous(C) and binary(B) outcomes."}			
	if (file.exists(paste0(train_file,".bim")) & file.exists(paste0(train_file,".bed")) & file.exists(paste0(train_file,".fam"))){} else {out1="The train file doesn't exist!"}
	if (file.exists(paste0(test_file,".bim")) & file.exists(paste0(test_file,".bed")) & file.exists(paste0(test_file,".fam"))){} else {out1="The test file doesn't exist!"}
	if (!file.exists(sum_stats_file)){out1="The summary statistic file does not exist!"} else {
		temp=fread(sum_stats_file,nrow=1)
		if (ncol(temp)==3){
			if (sum(colnames(temp) %in% c("V1","V2","V3"))==3){} else{
				if (sum(colnames(temp) %in% c("SNP","A1","Beta"))==3){} else {
					out1="The structure of sum_stats_file is wrong!"
				}
			} 
		} else {
			if (ncol(temp)>3){
				if (sum(colnames(temp) %in% c("SNP","A1","Beta"))<3){ 
					out1="The structure of sum_stats_file is wrong!"
				}
			} else {
				out1="The structure of sum_stats_file is wrong!"
			}
		} 
	}
	if (!LDblocks %in% c("EUR.hg19", "AFR.hg19", "ASN.hg19")) {out1="The LDblocks name is wrong!"}
	return(out1)
}


splitgenome2<-function (CHR, POS, ref.CHR, ref.breaks, details = T, right = TRUE)
{
    CHR <- as.character(CHR)
    ref.CHR <- as.character(ref.CHR)
    POS <- as.integer(POS)
    ref.breaks <- as.integer(ref.breaks)
    stopifnot(all(!is.na(POS)) && all(!is.na(ref.breaks)))
    stopifnot(all(POS >= 1) && all(ref.breaks >= 1))
    stopifnot(length(CHR) == length(POS))
    stopifnot(length(ref.CHR) == length(ref.breaks))
    chr <- (unique(CHR))
    chr.ref <- (unique(ref.CHR))
    included.chr <- chr %in% chr.ref
    if (!all(included.chr))
        stop("Some chromosomes were not defined in the reference. Make sure the notations are the same. e.g. 'chr1' vs 1 or chrX vs chr23.")
    levels <- character(0)
    results <- character(length(POS))
    Details <- data.frame()
    for (C in chr.ref) {
        breaks <- sort(unique(ref.breaks[ref.CHR == C]))
        if (breaks[1] > 1)
            breaks <- c(1, breaks)
        if (breaks[length(breaks)] < Inf)
            breaks <- c(breaks, Inf)
        cut <- cut(POS[CHR == C], include.lowest = T, breaks = breaks,
            right = right)
        levels <- c(levels, paste0(C, "_", levels(cut)))
        cut <- paste0(C, "_", as.character(cut))
        results[CHR == C] <- cut
        if (details) {
            df <- data.frame(chr = C, start = breaks[-length(breaks)],
                end = breaks[-1])
            Details <- rbind(Details, df)
        }
    }
    results <- factor(results, levels = levels)
    if (details) {
        Details$counts <- as.integer(table(results))
        attr(results, "details") <- Details
    }
   out.item=list(results,Details)
   return(out.item)
}

split_SNPandA1<-function(x){ 
	xtemp=strsplit(x,"_")[[1]]
	if (length(xtemp)==2){
		return(xtemp[1:2])
	} else {
		xn=length(xtemp);
		return(c(paste0(xtemp[1:(xn-1)],collapse="_"),xtemp[xn]))
	}
} 



PRStr_tuning<-function(Beta.all, ped, Covar_name,Y_name, Ytype, test_file){

    beta.info=Beta.all[,1:2]
    for (i in 9:ncol(Beta.all)){
      sdtemp=sd(Beta.all[,i],na.rm=T)
      if (sdtemp>1){
        Beta.all[,i:ncol(Beta.all)]=1;##print("fdsf")
      }
    }
    
    beta.all=Beta.all[,9:ncol(Beta.all)]/Beta.all$sd
    PRS.all=Calculate_PRS(test_file,beta.info,beta.all)
    
 
    out.item=data.frame("order"=NA,"R2"=NA)
    for (flag in 3:ncol(PRS.all)){
      out.item[flag-2,1]=flag-2
      resulttemp=PRS.all[,c(1,2,flag)]
      colnames(resulttemp)[3]="SCORESUM"
      if (Ytype=="C"){
        out.item[flag-2,2]=as.numeric(linear_result_generator(resulttemp,ped,Covar_name,Y_name))
      } else {
        out.item[flag-2,2]=as.numeric(logistic_result_generator(resulttemp,ped,Covar_name,Y_name))
      }
    }
    flag=which(out.item$R2==max(out.item$R2))[1]
    param_table=data.frame("lr"=c(0,rep(1/nrow(beta.all)*c(1,10,100,1000),each=15)),"iter"=c(0,rep(1:15,4)))   
    out.final=list()
    out.final$best.learning.rate=param_table$lr[flag]
    out.final$best.iteration=param_table$iter[flag]
    out.final$best.beta=cbind(beta.info,"beta"=beta.all[,flag]) 
    out.final$best.PRS=PRS.all[,c(1:3,flag+2)] ; colnames(out.final$best.PRS)[3:4]=c("PRS.NULL","PRS.TL")
    out.final$param_table=param_table 
    return(out.final)
}

