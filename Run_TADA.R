## Study: Rare Variant Analyses in Ancestrally Diverse Cohorts Reveal Novel ADHD Risk Genes
## Analysis: Run Transmission and De Novo Association Test (TADA)
## Author: Seulgi Jung

rm(list=ls())

library(epitools)
library(ggplot2)
library(webr)
library(dplyr)
library(ggsci)
library(wesanderson)
library(scales)
library(tidyverse)
library(ggpubr)
library(grid)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(magick)
library(ragg)
library(sjPlot)
library(svglite)
library(investr)

## read a input file including rare variant count per gene in each dataset
d <- read.table("input_file.txt", header=T)

## For TADA, Olfson + SPARK + Total Case-control
data <- data.frame(gene=d$gene, gene38=d$gene38, ENSG=d$ENSG, chromosome=d$chromosome, pLI=d$pLI, LOEUF=d$LOEUF, 
                    mut.ptv.g=d$mut.ptv.g, mut.misb.g=d$mut.misb.g, mut.misa.g=d$mut.misa.g, mut.syn.g=d$mut.syn.g, 
                    mut.ptv=d$mut.ptv, mut.misb=d$mut.misb, mut.misa=d$mut.misa, mut.syn=d$mut.syn, 
                    mut.mis.pathogenic=d$mut.mis.pathogenic, mut.mis.ambiguous=d$mut.mis.ambiguous,
                    mean_am_pathogenicity=d$mean_am_pathogenicity, avg_alphamis_pathogenic=d$avg_alphamis_pathogenic,
                    avg_alphamis_ambiguous=d$avg_alphamis_ambiguous,
                    dn.ptv = d$spark_ptv_prob + d$olfson_ptv_prob,
                    dn.misb = d$spark_misb_prob + d$olfson_misb_prob,  
                    dn.misa = d$spark_misa_prob + d$olfson_misa_prob,
                    dn.mis = d$spark_mis_prob + d$olfson_mis_prob,
                    dn.mis_total = d$spark_misb_prob + d$olfson_misb_prob + d$spark_misa_prob + d$olfson_misa_prob + d$spark_mis_prob + d$olfson_mis_prob,
                    dn.mis.pathogenic = d$spark_mis_pathogenic_prob + d$olfson_mis_pathogenic_prob,  
                    dn.mis.ambiguous = d$spark_mis_ambiguous_prob + d$olfson_mis_ambiguous_prob,
                    dn.mis.benign = d$spark_mis_benign_prob + d$olfson_mis_benign_prob,
                    dn.syn = d$spark_syn_prob + d$olfson_syn_prob,
                    
                    dn.ptv.sib = d$spark_ptv_sib + d$olfson_ptv_sib, 
                    dn.misb.sib = d$spark_misb_sib + d$olfson_misb_sib, 
                    dn.misa.sib = d$spark_misa_sib + d$olfson_misa_sib,
                    dn.mis.sib = d$spark_mis_sib + d$olfson_mis_sib,
                    dn.mis_total.sib = d$spark_misb_sib + d$olfson_misb_sib + d$spark_misa_sib + d$olfson_misa_sib + d$spark_mis_sib + d$olfson_mis_sib,
                    dn.mis.pathogenic.sib = d$spark_mis_pathogenic_sib + d$olfson_mis_pathogenic_sib, 
                    dn.mis.ambiguous.sib = d$spark_mis_ambiguous_sib + d$olfson_mis_ambiguous_sib,
                    dn.mis.benign.sib = d$spark_mis_benign_sib + d$olfson_mis_benign_sib,
                    dn.syn.sib = d$spark_syn_sib + d$olfson_syn_sib,
                    
                    cc.ptv.case = d$allofus_ptv_case + d$kyle_ptv_case,
                    cc.misb.case = d$allofus_misb_case + d$kyle_misb_case,
                    cc.misa.case = d$allofus_misa_case,
                    cc.mis.case = d$allofus_mis_case,
                    cc.mis.pathogenic.case = d$allofus_mis_pathogenic_case,
                    cc.mis.ambiguous.case = d$allofus_mis_ambiguous_case,
                    cc.mis.benign.case = d$allofus_mis_benign_case,
                    cc.syn.case = d$allofus_syn_case + d$kyle_syn_case,
                    cc.syn.case.allofus = d$allofus_syn_case,
                    
                    cc.ptv.control = d$allofus_ptv_control + d$kyle_ptv_control,
                    cc.misb.control = d$allofus_misb_control + d$kyle_misb_control,
                    cc.misa.control = d$allofus_misa_control,
                    cc.mis.control = d$allofus_mis_control,
                    cc.mis.pathogenic.control = d$allofus_mis_pathogenic_control,
                    cc.mis.ambiguous.control = d$allofus_mis_ambiguous_control,
                    cc.mis.benign.control = d$allofus_mis_benign_control,
                    cc.syn.control = d$allofus_syn_control + d$kyle_syn_control,
                    cc.syn.control.allofus = d$allofus_syn_control
                    )



## Number of samples 
## trio datasets
spark_prob = 902
spark_sib = 3508

olfson_prob = 147
olfson_sib = 780

## case-control datasets
allofus_case = 4856
allofus_control = 24108

kyle_case = 3477
kyle_control = 5002

## total
n_prob <- olfson_prob + spark_prob
n_sib <- olfson_sib + spark_sib

n1_case = allofus_case+kyle_case
n1_control = allofus_control + kyle_control

n2_case = allofus_case
n2_control = allofus_control 



## Functions
evidence.null.cc <- function(x.cc, n.cc, rho0, nu0) {
  marglik0.ctrl.log <- log(dnbinom(x.cc$cn, rho0, nu0/(nu0+n.cc$cn)))
  marglik0.case.log <- log(dnbinom(x.cc$ca, rho0+x.cc$cn, (nu0+n.cc$cn)/(nu0+n.cc$cn+n.cc$ca)))
  marglik0.log <- marglik0.ctrl.log + marglik0.case.log
  return (list(cn=exp(marglik0.ctrl.log), ca=exp(marglik0.case.log), total=exp(marglik0.log)))
}

evidence.alt.cc <- function(x.cc, n.cc, gamma.cc, beta.cc, rho1, nu1, q.lower=0, q.upper=Inf) {
  integrand <- function(q) {
    return (dnbinom(x.cc$ca, gamma.cc*beta.cc, beta.cc) * dgamma(q, rho1+x.cc$cn, nu1+n.cc$cn))
  }
  
  marglik1.ctrl <- dnbinom(x.cc$cn, rho1, nu1/(nu1+n.cc$cn))
  marglik1.case <- dnbinom(x.cc$ca, rho1+x.cc$cn, (nu1+n.cc$cn)/(nu1+n.cc$ca*gamma.cc+n.cc$cn))
  marglik1 <- marglik1.ctrl * marglik1.case
  return (list(cn=marglik1.ctrl, ca=marglik1.case, total=marglik1))
}

table.BF.dn.wrapper <- function(i.x,x,n.dn,mu,gamma.dn,beta.dn){
  BF=bayes.factor.dn(x[i.x],n.dn=n.dn,mu=mu,gamma.dn=gamma.dn,beta.dn=beta.dn)
  return(BF)
}

bayes.factor.dn <- function(x.dn,n.dn,mu,gamma.dn,beta.dn){
  marg.lik0 <- dpois(x.dn, 2*n.dn*mu)
  marg.lik1 <- dnbinom(x.dn, gamma.dn*beta.dn, beta.dn/(beta.dn+2*n.dn*mu))
  BF <- marg.lik1/marg.lik0
  return (BF=BF)
}

bayes.factor.cc <- function(x.cc, n.cc, gamma.cc, beta.cc,rho1,nu1,rho0,nu0){
  marglik0.cc <- evidence.null.cc(x.cc, n.cc, rho0, nu0)
  marglik1.cc <- evidence.alt.cc(x.cc, n.cc, gamma.cc, beta.cc, rho1, nu1)
  BF.cn <- marglik1.cc$cn / marglik0.cc$cn
  BF.ca <- marglik1.cc$ca / marglik0.cc$ca
  BF <- BF.cn * BF.ca
  return(BF=BF)
}



## Calculate Bayes factors (BFs): de novo
gamma_olfson <- read.table("dn_gamma_pi0.01_Olfson.txt", header=T)

tada$gamma_dn_ptv = gamma_olfson$gamma_dn_ptv
tada$gamma_dn_misb <- gamma_olfson$gamma_dn_misb_fix
tada$gamma_dn_misa <- gamma_olfson$gamma_dn_misa_fix

cat0 = c('ptv', 'misa', 'misb')
gamma_dn = matrix(NA, nrow(tada), length(cat0))
colnames(gamma_dn) = cat0
gamma_dn[,'ptv'] = tada$gamma_dn_ptv
gamma_dn[,'misb'] = tada$gamma_dn_misb
gamma_dn[,'misa'] = tada$gamma_dn_misa


beta_dn = 0.2

BF_dn = matrix(NA, nrow(tada), length(cat0))
colnames(BF_dn) = cat0
for(j in (cat0)) {
  # print(j)
  BF_dn[,j] = sapply(1:nrow(tada), function(i) bayes.factor.dn(x.dn = tada[i, paste0('dn.', j)], n.dn = n_prob, mu = tada[i, paste0('mut.', j)], 
                                                               beta.dn = beta_dn, gamma.dn = gamma_dn[i,j]))
}
head(BF_dn)

BF_dn_ptv <- BF_dn[,1]
BF_dn_misb <- BF_dn[,3]
BF_dn_misa <- BF_dn[,2]


## Calculating BFs: case-control
gamma_cc_allofus <- read.table("cc_gamma_pi0.01_allofus.txt", header=T)

tada$gamma_cc_ptv = gamma_cc_allofus$gamma_cc_ptv
tada$gamma_cc_misb <- gamma_cc_allofus$gamma_cc_misb_fix
tada$gamma_cc_misa <- gamma_cc_allofus$gamma_cc_misa_fix

tada$rho.ptv <- gamma_cc_allofus$rho.ptv
tada$rho.misb <- gamma_cc_allofus$rho.misb
tada$rho.misa <- gamma_cc_allofus$rho.misa


##  ccPTV
n.cc = data.frame(ca=n1_case, cn=n1_control)
BF_cc_ptv = 
  sapply((1:nrow(tada)), 
         function(i) bayes.factor.cc(x.cc= data.frame(ca=tada[i, paste0('cc.ptv.case')], cn=tada[i, paste0('cc.ptv.control')]), 
                                     n.cc = n.cc, gamma.cc=tada$gamma_cc_ptv[i], beta.cc=NULL, rho1=tada$rho.ptv[i], 
                                     nu1=nu, rho0=tada$rho.ptv[i], nu0=nu))


## ccMisB
n.cc = data.frame(ca=n1_case, cn=n1_control)

BF_cc_misb = 
  sapply((1:nrow(tada)), 
         function(i) bayes.factor.cc(x.cc= data.frame(ca=tada[i, paste0('cc.misb.case')], cn=tada[i, paste0('cc.misb.control')]), 
                                     n.cc = n.cc, gamma.cc=tada$gamma_cc_misb[i], beta.cc=NULL, rho1=tada$rho.misb[i], 
                                     nu1=nu, rho0=tada$rho.misb[i], nu0=nu))

## ccMisA 
n.cc2 = data.frame(ca=n2_case, cn=n2_control)
BF_cc_misa = sapply((1:nrow(tada)), function(i) bayes.factor.cc(x.cc= data.frame(ca=tada[i, paste0('cc.misa.case')], cn=tada[i, paste0('cc.misa.control')]), 
                                                                n.cc = n.cc2, gamma.cc=tada$gamma_cc_misa[i], beta.cc=NULL, rho1=tada$rho.misa[i], 
                                                                nu1=nu, rho0=tada$rho.misa[i], nu0=nu))



## Multiply BFs
BF_dn_ptv[BF_dn_ptv<1] <- 1
BF_dn_misb[BF_dn_misb<1] <- 1
BF_dn_misa[BF_dn_misa<1] <- 1
BF_cc_ptv[BF_cc_ptv<1] <- 1
BF_cc_misb[BF_cc_misb<1] <- 1
BF_cc_misa[BF_cc_misa<1] <- 1

tada$BF_dn_ptv <- BF_dn_ptv
tada$BF_dn_misb <- BF_dn_misb
tada$BF_dn_misa <- BF_dn_misa
tada$BF_cc_ptv <- BF_cc_ptv
tada$BF_cc_misb <- BF_cc_misb
tada$BF_cc_misa <- BF_cc_misa


BF <- cbind(pmax(BF_dn_ptv*BF_cc_ptv, 1),
                pmax(BF_dn_misb*BF_cc_misb, 1),
                pmax(BF_dn_misa*BF_cc_misa, 1))

tada$BF = BF

bf <- data.frame(BF)
head(bf)


## calculation of false discovery rate (FDR)
Bayesian.FDR <- function(BF, pi0) {
  i.order=order(BF, decreasing = T)
  BF=BF[i.order]
 
  pi <- 1-pi0
  q <- pi*BF/(1-pi+pi*BF)
  q0 <- 1 - q 
  
  FDR=cumsum(q0)/(1:length(BF))
  FDR[i.order]=FDR
  
  return (FDR=FDR)
}


## calculation of P value
q2p <- function(x){
  ec <- x
  ec <- sort(ec)
  ec <- ec*(1:length(ec))/(length(ec))/max(ec)    
  return(p=ec[names(x)])
}


## Add FDR and P value to gene
tada$qval = Bayesian.FDR(apply(BF, 1, prod), pi0 = 1 - 0.01)
fdr <- tada$qval
tada$qval_final <- ifelse(is.na(tada$qval), 1, tada$qval)

qval <- Bayesian.FDR(apply(BF, 1, prod), pi0 = 1 - 0.01)
names(qval) <- tada$ENSG
tada$pval <- q2p(qval)


write.table(tada, "TADA_result.txt", col.names=T, row.names=F, quote=F, sep="\t")

