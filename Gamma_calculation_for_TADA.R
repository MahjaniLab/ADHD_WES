## Study: Rare Variant Analyses in Ancestrally Diverse Cohorts Reveal Novel ADHD Risk Genes
## Analysis: Calculating gamma to perform the Transmission and De Novo Association Test (TADA)
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

## For calculating gamma
data <- data.frame(gene=d$gene, gene38=d$gene38, ENSG=d$ENSG, chromosome=d$chromosome, pLI=d$pLI, LOEUF=d$LOEUF, 
                     mut.ptv.g=d$mut.ptv.g, mut.misb.g=d$mut.misb.g, mut.misa.g=d$mut.misa.g, mut.syn.g=d$mut.syn.g, 
                     mut.ptv=d$mut.ptv, mut.misb=d$mut.misb, mut.misa=d$mut.misa, mut.syn=d$mut.syn, 
                     mut.mis.pathogenic=d$mut.mis.pathogenic, mut.mis.ambiguous=d$mut.mis.ambiguous,
                     mean_am_pathogenicity=d$mean_am_pathogenicity, avg_alphamis_pathogenic=d$avg_alphamis_pathogenic,
                     avg_alphamis_ambiguous=d$avg_alphamis_ambiguous,
                     dn.ptv = d$olfson_ptv_prob,
                     dn.misb = d$olfson_misb_prob,
                     dn.misa = d$olfson_misa_prob,
                     dn.mis = d$olfson_mis_prob,
                     dn.syn = d$olfson_syn_prob,
                     
                     dn.ptv.sib = d$olfson_ptv_sib, 
                     dn.misb.sib = d$olfson_misb_sib,
                     dn.misa.sib = d$olfson_misa_sib,
                     dn.mis.sib = d$olfson_mis_sib,
                     dn.syn.sib = d$olfson_syn_sib,
                     
                     
                     cc.ptv.case = d$allofus_ptv_case ,
                     cc.misb.case = d$allofus_misb_case ,
                     cc.misa.case = d$allofus_misa_case,
                     cc.mis.case = d$allofus_mis_case,
                     cc.mis.pathogenic.case = d$allofus_mis_pathogenic_case,
                     cc.mis.ambiguous.case = d$allofus_mis_ambiguous_case,
                     cc.mis.benign.case = d$allofus_mis_benign_case,
                     cc.syn.case = d$allofus_syn_case,
                     
                     cc.ptv.control = d$allofus_ptv_control ,
                     cc.misb.control = d$allofus_misb_control ,
                     cc.misa.control = d$allofus_misa_control,
                     cc.mis.control = d$allofus_mis_control,
                     cc.mis.pathogenic.control = d$allofus_mis_pathogenic_control,
                     cc.mis.ambiguous.control = d$allofus_mis_ambiguous_control,
                     cc.mis.benign.control = d$allofus_mis_benign_control,
                     cc.syn.control = d$allofus_syn_control
                     
)



## sample number of trio datset
olfson_prob = 147
olfson_sib = 780
sample.n <- c(olfson_prob, olfson_sib)

## Adjustment of mutation rate based on the data
actual <- apply(data[,c("dn.ptv.sib", "dn.misb.sib", "dn.misa.sib")], 2, sum)
expected <- round(apply(data[,c("mut.ptv.g", "mut.misb.g", "mut.misa.g")]*(sample.n[2]), 2, sum)*2)
cal <- actual/expected

data$mut.ptv.g <- data$mut.ptv.g*cal[1]
data$mut.misb.g <- data$mut.misb.g*cal[2]
data$mut.misa.g <- data$mut.misa.g*cal[3]

### Over-write mutation rates with modified gnomad rates
data$mut.ptv <- data$mut.ptv.g
data$mut.misb <- data$mut.misb.g
data$mut.misa <- data$mut.misa.g

data[is.na(data)] <- 0
tada <- data

####################################################

## sample number of case-control dataset
allofus_case = 4856
allofus_control = 24108

n_prob <- olfson_prob
n_sib <- olfson_sib

n1_case = allofus_case
n1_control = allofus_control
n2_case = allofus_case
n2_control = allofus_control 



# gamma calculation of de novo protein-truncating variant (ptv)

tada$loeuf <- tada$LOEUF/2

n.trio = data.frame(ca=n_prob, cn=n_sib)

bin = 5
pi = rep(0.01, bin) 
nbins = bin
loeuf_max = 1
nu = 1e2
bins = 1:nbins/nbins
uni = (quantile(tada$loeuf, bins - bins[1]/2))
uni*2

interval = cbind(quantile(tada$loeuf, bins - bins[1]), quantile(tada$loeuf, bins))
interval

lambda_dn_ptv = ngene = rep(NA, length(uni))
for(k in 1:length(uni)) {
  input <- tada[which(tada$loeuf >= interval[k,1] & tada$loeuf <= interval[k,2]),]
  ngene[k] <- dim(input)[1]
  case_ptv <- sum(input$dn.ptv)
  control_ptv <- sum(input$dn.ptv.sib)
  case_syn <- sum(input$dn.syn)
  control_syn <- sum(input$dn.syn.sib)
  lambda_dn_ptv[k] <- (case_ptv/case_syn)/(control_ptv/control_syn)
  print(lambda_dn_ptv[k])
}

lambda_dn_ptv
lambda_dn_ptv[4:5] <- 1

y_lambda = pmax(1, lambda_dn_ptv)

x = qnorm(uni)
x

logisticModelSS_lambda <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = list(y = y_lambda, x = x))
pred_gamma_pli_lambda = predFit(logisticModelSS_lambda, interval="prediction", newdata = data.frame(x = qnorm(tada$loeuf)))
lambda_dn_ptv = pmax(1, pmax(1, pred_gamma_pli_lambda[, "fit"]))
tada$lambda_dn_ptv = lambda_dn_ptv

plot(y=tada$lambda_dn_ptv, x=tada$LOEUF[!is.na(tada$LOEUF)], ylim=c(0.5, 5), ylab="Lambda - de novo PTV", xlab="LOEUF")
points(y=y_lambda, x=uni*2, col=2)

par(new = TRUE)
data <- data.frame(x = uni * 2, y = y_lambda)
barplot(height = data$y, 
        names.arg = data$x[!is.na(data$x)], 
        ylim = c(0, 5), 
        axes = FALSE,  
        col = rgb(1, 0, 0, 0.5),  
        border = NA)  

gamma_dn_ptv <- 1 + (lambda_dn_ptv - 1) / pi
gamma_dn_ptv

## smmothing gamma
y = pmax(1, gamma_dn_ptv)
y

x = qnorm(uni)
x

plot(x, y, main = "Data", xlab = "x", ylab = "y")

logisticModelSS <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = list(y = y, x = x))
pred_gamma_pli = predFit(logisticModelSS, interval="prediction", newdata = data.frame(x = qnorm(tada$loeuf)))
gamma_dn_ptv = pmax(1, pmax(1, pred_gamma_pli[, "fit"]))
tada$gamma_dn_ptv = gamma_dn_ptv

plot(y=tada$gamma_dn_ptv, x=tada$LOEUF[!is.na(tada$LOEUF)], ylim=c(0.5, 400), ylab="Gamma de novo PTV", xlab="LOEUF")


## smmothing gamma de novo misB
tada$alphamis <- 1-tada$mean_am_pathogenicity

bin = 5
pi = rep(0.01, bin) ## 
nbins = bin
loeuf_max = 1
nu = 1e2

bins = 1:nbins/nbins

uni = (quantile(tada$alphamis, bins - bins[1]/2))
uni

interval = cbind(quantile(tada$alphamis, bins - bins[1]), quantile(tada$alphamis, bins))
interval

lambda_dn_misB = ngene = rep(NA, length(uni))
for(k in 1:length(uni)) {
  input <- tada[which(tada$alphamis >= interval[k,1] & tada$alphamis <= interval[k,2]),]
  ngene[k] <- dim(input)[1]
  case_ptv <- sum(input$dn.misb)
  control_ptv <- sum(input$dn.misb.sib)
  case_syn <- sum(input$dn.syn)
  control_syn <- sum(input$dn.syn.sib)
  lambda_dn_misB[k] <- (case_ptv/case_syn)/(control_ptv/control_syn)
  print(lambda_dn_misB[k])
  #print(case_ptv/case_syn)
}

lambda_dn_misB
lambda_dn_misB[3:5] <- 1

gamma_dn_misB <- 1 + (lambda_dn_misB - 1) / pi
gamma_dn_misB

y = pmax(1, gamma_dn_misB)
y

x = qnorm(uni)
x

plot(x, y, main = "Data", xlab = "x", ylab = "y")

logisticModelSS <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = list(y = y, x = x))
pred_gamma_pli = predFit(logisticModelSS, interval="prediction", newdata = data.frame(x = qnorm(tada$alphamis)))
gamma_dn_misB = pmax(1, pmax(1, pred_gamma_pli[, "fit"]))
tada$gamma_dn_misb = gamma_dn_misB

plot(y=tada$gamma_dn_misb, x=tada$alphamis[!is.na(tada$alphamis)], ylim=c(0.5, 400), ylab="Case-control gamma MisB", xlab="1-alphamis")


## smmothing gamma de novo misA
tada$alphamis <- 1-tada$mean_am_pathogenicity

bin = 4
pi = rep(0.01, bin) ## 
nbins = bin
loeuf_max = 1
nu = 1e2

bins = 1:nbins/nbins

uni = (quantile(tada$alphamis, bins - bins[1]/2))
uni

interval = cbind(quantile(tada$alphamis, bins - bins[1]), quantile(tada$alphamis, bins))
interval

lambda_dn_misA = ngene = rep(NA, length(uni))
for(k in 1:length(uni)) {
  input <- tada[which(tada$alphamis >= interval[k,1] & tada$alphamis <= interval[k,2]),]
  ngene[k] <- dim(input)[1]
  case_ptv <- sum(input$dn.misa)
  control_ptv <- sum(input$dn.misa.sib)
  case_syn <- sum(input$dn.syn)
  control_syn <- sum(input$dn.syn.sib)
  lambda_dn_misA[k] <- (case_ptv/case_syn)/(control_ptv/control_syn)
  print(lambda_dn_misA[k])
  #print(case_ptv/case_syn)
}

lambda_dn_misA
lambda_dn_misA[3:4] <- 1

gamma_dn_misA <- 1 + (lambda_dn_misA - 1) / pi
gamma_dn_misA

y = pmax(1, gamma_dn_misA)
y

x = qnorm(uni)
x

plot(x, y, main = "Data", xlab = "x", ylab = "y")

logisticModelSS <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = list(y = y, x = x))
pred_gamma_pli = predFit(logisticModelSS, interval="prediction", newdata = data.frame(x = qnorm(tada$alphamis)))
gamma_dn_misA = pmax(1, pmax(1, pred_gamma_pli[, "fit"]))
tada$gamma_dn_misa <- gamma_dn_misA
plot(y=tada$gamma_dn_misa, x=tada$alphamis[!is.na(tada$alphamis)], ylim=c(0.5, 100), ylab="Case-control gamma MisB", xlab="1-alphamis")

################################################################################


## Save the gamma of de novo variants
gamma = c(gene = NA, gamma_dn_ptv = NA, gamma_dn_misB = NA, gamma_dn_misA = NA, gamma_dn_misb_fix = NA, gamma_dn_misa_fix = NA) 
gamma['gene'] <- list(tada$gene)
gamma['gamma_dn_ptv'] <- list(gamma_dn_ptv)
gamma['gamma_dn_misB'] <- list(gamma_dn_misB)
gamma['gamma_dn_misA'] <- list(gamma_dn_misA)
gamma_olfson <- data.frame(gamma)
head(gamma_olfson)
write.table(gamma_olfson, "dn_gamma_pi0.01_Olfson.txt", col.names=T, row.names=F, quote=F, sep="\t")



######### Case-control #########################################################
######### gamma ptv

tada$loeuf <- tada$LOEUF/2

bin = 120
pi = rep(0.01, bin)
nbins = bin
loeuf_max = 1
nu = 1e2
n.cc = data.frame(ca=n1_case, cn=n1_control)

bins = 1:nbins/nbins
uni = (quantile(tada$loeuf, bins - bins[1]/2))
uni

interval = cbind(quantile(tada$loeuf, bins - bins[1]), quantile(tada$loeuf, bins))
interval

lambda0 = ngene = rep(NA, length(uni))
for(k in 1:length(uni)) {
  input <- tada[which(tada$loeuf >= interval[k,1] & tada$loeuf <= interval[k,2]),]
  ngene[k] <- dim(input)[1]
  case_ptv <- sum(input$cc.ptv.case)
  control_ptv <- sum(input$cc.ptv.control)
  case_syn <- sum(input$cc.syn.case)
  control_syn <- sum(input$cc.syn.control)
  lambda0[k] <- (case_ptv/case_syn)/(control_ptv/control_syn)
  print(lambda0[k])
  }

lambda0
lambda0[30:120] <- 1

gamma0 <- 1 + (lambda0 - 1) / pi
gamma0

## smoothing gamma
y = pmax(1, gamma0)
y

x = qnorm(uni)
x

plot(x, y, main = "Data", xlab = "x", ylab = "y")

logisticModelSS <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = list(y = y, x = x))
pred_gamma_pli = predFit(logisticModelSS, interval="prediction", newdata = data.frame(x = qnorm(tada$loeuf)))
gamma_cc_ptv = pmax(1, pmax(1, pred_gamma_pli[, "fit"]))
tada$gamma_cc_ptv <- gamma_cc_ptv
plot(y=tada$gamma_cc_ptv, x=tada$LOEUF[!is.na(tada$LOEUF)], ylim=c(0.5, 200), ylab="Case-control gamma PTV", xlab="LOEUF")

rho.ptv = nu * sum(tada$cc.ptv.control) / (2*nrow(tada)*(n.cc[2]))
tada$rho.ptv = as.numeric(rho.ptv)*tada$mut.ptv.g/mean(tada$mut.ptv.g)



## gamma misB
tada$alphamis <- 1-tada$mean_am_pathogenicity

bin = 6
pi = rep(0.01, bin)
nbins = bin
loeuf_max = 1
nu = 1e2
n.cc = data.frame(ca=n1_case, cn=n1_control)

bins = 1:nbins/nbins

tada$cc_syn = (tada$cc.syn.case >0) | (tada$cc.syn.control >0)

uni = (quantile(tada$loeuf[tada$cc_syn], bins - bins[1]/2))
uni

interval = cbind(quantile(tada$loeuf[tada$cc_syn], bins - bins[1]), quantile(tada$loeuf[tada$cc_syn], bins))
interval

lambda_misB = ngene = rep(NA, length(uni))
for(k in 1:length(uni)) {
  input <- tada[which(tada$alphamis >= interval[k,1] & tada$alphamis <= interval[k,2]),]
  ngene[k] <- dim(input)[1]
  case_misb <- sum(input$cc.misb.case)
  control_misb <- sum(input$cc.misb.control)
  case_syn <- sum(input$cc.syn.case)
  control_syn <- sum(input$cc.syn.control)
  lambda_misB[k] <- (case_misb/case_syn)/(control_misb/control_syn)
  print(lambda_misB[k])
  #print(case_ptv/case_syn)
}

lambda_misB
lambda_misB[5:6] <- 1

gamma_misB <- 1 + (lambda_misB - 1) / pi
gamma_misB

## smoothing gamma
y = pmax(1, gamma_misB)
y

x = qnorm(uni)
x

plot(x, y, main = "Data", xlab = "x", ylab = "y")

logisticModelSS <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = list(y = y, x = x))
pred_gamma_pli = predFit(logisticModelSS, interval="prediction", newdata = data.frame(x = qnorm(tada$alphamis)))
gamma_cc_misB = pmax(1, pmax(1, pred_gamma_pli[, "fit"]))
tada$gamma_cc_misb <- gamma_cc_misB
plot(y=tada$gamma_cc_misb, x=tada$alphamis[!is.na(tada$alphamis)], ylim=c(0.5, 50), ylab="Case-control gamma MisB", xlab="1-alphamis")

rho.misb = nu * sum(tada$cc.misb.control) / (2*nrow(tada)*(n.cc[2]))
tada$rho.misb = as.numeric(rho.misb)*tada$mut.misb.g/mean(tada$mut.misb.g)



## gamma misA
tada$alphamis <- 1-tada$mean_am_pathogenicity

bin = 6
pi = rep(0.01, bin)
nbins = bin
loeuf_max = 1
nu = 1e2
n.cc = data.frame(ca=n2_case, cn=n2_control)

bins = 1:nbins/nbins

uni = (quantile(tada$alphamis, bins - bins[1]/2))
uni[uni == 1] <- 0.8
uni

interval = cbind(quantile(tada$alphamis, bins - bins[1]), quantile(tada$alphamis, bins))
interval

lambda_misA = ngene = rep(NA, length(uni))
for(k in 1:length(uni)) {
  input <- tada[which(tada$alphamis >= interval[k,1] & tada$alphamis <= interval[k,2]),]
  ngene[k] <- dim(input)[1]
  case_misa <- sum(input$cc.misa.case)
  control_misa <- sum(input$cc.misa.control)
  case_syn <- sum(input$cc.syn.case)
  control_syn <- sum(input$cc.syn.control)
  lambda_misA[k] <- (case_misa/case_syn)/(control_misa/control_syn)
  print(lambda_misA[k])
  }

lambda_misA
lambda_misA[5:6] <- 1

gamma_misA <- 1 + (lambda_misA - 1) / pi
gamma_misA

## Smoothing gamma
y = pmax(1, gamma_misA)
y

x = qnorm(uni)
x

plot(x, y, main = "Data", xlab = "x", ylab = "y")

logisticModelSS <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = list(y = y, x = x))
pred_gamma_pli = predFit(logisticModelSS, interval="prediction", newdata = data.frame(x = qnorm(tada$alphamis)))
gamma_cc_misA = pmax(1, pmax(1, pred_gamma_pli[, "fit"]))
tada$gamma_cc_misA = gamma_cc_misA

plot(y=tada$gamma_cc_misA, x=tada$alphamis[!is.na(tada$alphamis)], ylim=c(0.5, 50), ylab="Case-control gamma MisB", xlab="1-alphamis")


n.cc2 = data.frame(ca=n2_case, cn=n2_control)
rho.misa = nu * sum(tada$cc.misa.control) / (2*nrow(tada)*(n.cc2[2]))
tada$rho.misa = as.numeric(rho.misa)*tada$mut.misa.g/mean(tada$mut.misa.g)



################################################################################
## Save the case-control gamma
gamma_cc = c(gene = NA, gamma_cc_ptv = NA, gamma_cc_misB = NA, gamma_cc_misA = NA, gamma_cc_misb_fix = NA, gamma_cc_misa_fix = NA,
             rho.ptv = NA, rho.misb=NA, rho.misa=NA) 

gamma_cc['gene'] <- list(tada$gene)
gamma_cc['gamma_cc_ptv'] <- list(gamma_cc_ptv)
gamma_cc['gamma_cc_misB'] <- list(gamma_cc_misB)
gamma_cc['gamma_cc_misA'] <- list(gamma_cc_misA)
gamma_cc['rho.ptv'] <- rho.ptv
gamma_cc['rho.misb'] <- rho.misb
gamma_cc['rho.misa'] <- rho.misa
gamma_cc_allofus <- data.frame(gamma_cc)
head(gamma_cc_allofus)
write.table(gamma_cc_allofus, "cc_gamma_pi0.01_allofus.txt", col.names=T, row.names=F, quote=F, sep="\t")

