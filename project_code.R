
# Printing pdf figures
printfig <- function(name, height = 6, ...){
  path <- paste('~/stk4080//latex/figures/', name, '.pdf', sep = '')
  pdf(file = path, height = height, ...)
  par(mai = c(0.6, 0.6, 0.1, 0.1)) # manipulate inner margins of subplots
  par(mgp = c(1.8, 0.7, 0)) # manipulate spacing of axis ticks, labels, text
  par(omi = c(0, 0, 0, 0)) # manipulate outer margin of full plot
}

pbc <- read.csv("~/stk4080/pbc3A.txt", sep="")
head(pbc)
# sapply(pbc, class)
# pbc$dod = as.factor(pbc$dod)
# pbc$beh = as.factor(pbc$beh)
# pbc$kjonn = as.factor(pbc$kjonn)
# plot(pbc$beh)
# plot(pbc$beh, pbc$dod)
# plot(pbc)
# plot(pbc$kjonn, pbc$ald)
?cox.zph 

# Not unique death / censor times
# We have tied data because of rounding to nearest day.
length(pbc$dager) 
length(unique(pbc$dager))
sort(pbc$dager)
?coxph
#----------
# a) Enkle univariate analyser:
# Først utfører du enkle univariate analyser ved hjelp av Kaplan-Meier plot og 
# logrank tester hvor effekten av de fem kovariatene studeres hver for seg. 
# For de numeriske kovariatene alder, bilirubin og albumin ma du selv velge en 
# passende gruppering.

#---------------
# Kaplan-Meier
#---------------
library(survival)
printfig('km_beh', height = 5, width=5)
layout(matrix(c(rep(1, 3), rep(2, 2)), nrow = 5, ncol = 1))
fit = survfit(Surv(dager, dod)~beh, data = pbc)
survfit(Surv(pbc$dager, pbc$dod)~pbc$beh)
plot(fit, lty=1:2, ylab='survival', xlab='days')
legend(1,0.2,c("treatment","placebo"), lty=1:2, bty='n')
grid()
plot(as.factor(pbc$beh), main='', ylab='frequency', xlab='treatment', col='white',
     names.arg=c('treatment', 'placebo'))
dev.off()

survdiff(Surv(dager, dod)~beh, data = pbc)

library(dplyr)
df = tbl_df(pbc)
glimpse(df)
df %>% group_by(kjonn) %>% summarise(avg = mean(dager))

str_lab = function(gr){
  gr = as.character(gr)
  lab = c(paste('-', gr[2]))
  for (i in 2:(length(gr)-1)){
    lab[i] = paste(gr[i], '-', gr[i+1])
  }
  lab[i+1] = paste(gr[i+1], '-')
  lab
}

#----
printfig('km_ald', height = 5, width=5)
layout(matrix(c(rep(1, 3), rep(2, 2)), nrow = 5, ncol = 1))
pbc$ald_groups = 0
gr = c(0, 40, 50, 60)
for (i in gr) {
  pbc$ald_groups[pbc$ald > i] = i
}
fit = survfit(Surv(dager, dod)~ald_groups, data = pbc)
plot(fit, lty=1:length(gr), ylab='survival', xlab='days')
legend(1,0.3, str_lab(gr), lty=1:length(gr), bty='n')
grid()
hist(pbc$ald, breaks = 20, main='', ylab='frequency', xlab='age')
dev.off()

survdiff(Surv(dager, dod)~beh, data = pbc)
survdiff(Surv(dager, dod)~ald_groups, data = pbc)


#----
printfig('km_kjonn', height = 5, width=5)
layout(matrix(c(rep(1, 3), rep(2, 2)), nrow = 5, ncol = 1))
fit = survfit(Surv(dager, dod)~kjonn, data = pbc)
plot(fit, lty=1:2, ylab='survival', xlab='days')
legend(1,0.2,c("male","female"), lty=1:2, bty='n')
grid()
plot(as.factor(pbc$kjonn), main='', ylab='frequency', xlab='gender', col='white',
     names.arg=c('male', 'female'))
dev.off()

survdiff(Surv(dager, dod)~kjonn, data = pbc)

#----
printfig('km_bil', height = 5, width=5)
layout(matrix(c(rep(1, 3), rep(2, 2)), nrow = 5, ncol = 1))
pbc$bil_groups = 0
gr = c(0, 1, 3, 7, 15)
for (i in gr) {
  pbc$bil_groups[pbc$bil > i] = i
}
fit = survfit(Surv(dager, dod)~bil_groups, data = pbc)
plot(fit, lty=1:length(gr), ylab='survival', xlab='days')
legend(1,0.4, str_lab(gr), lty=1:length(gr), bty='n')
grid()
hist(pbc$bil, breaks = 20, main='', ylab='frequency', xlab='bilirubin (mg/dl)')
dev.off()

survdiff(Surv(dager, dod)~bil_groups, data = pbc)

#----
printfig('km_alb', height = 5, width=5)
# par(mfrow = c(1, 1))
# par(mfrow = c(2, 1), mai = c(1, 1, 0.1, 0.1))
layout(matrix(c(rep(1, 3), rep(2, 2)), nrow = 5, ncol = 1))
pbc$alb_groups = 0
gr = c(0, 3, 3.5, 4)
for (i in gr) {
  pbc$alb_groups[pbc$alb > i] = i
}
fit = survfit(Surv(dager, dod)~alb_groups, data = pbc)
plot(fit, lty=1:length(gr), ylab='survival', xlab='days')
legend(1, 0.4, str_lab(gr), lty=1:length(gr), bty='n')
grid()

hist(pbc$alb, breaks = 20, main = '', ylab = 'frequency', xlab = 'albumin (mg/dl)')
dev.off()

# Log-rank tests
survdiff(Surv(dager, dod)~beh, data = pbc)
survdiff(Surv(dager, dod)~ald_groups, data = pbc)
survdiff(Surv(dager, dod)~kjonn, data = pbc)
survdiff(Surv(dager, dod)~bil_groups, data = pbc)
survdiff(Surv(dager, dod)~alb_groups, data = pbc)

#############################################################
# b) Univariat regresjon:
# Deretter utfører du univariate Cox-regresjonsanalyser for hver av kovariatene. 
# I denne sammenhengen ma du spesielt vurdere hvordan de numeriske kovariatene
# bør kodes.

#-----------
# Treatment
cox = coxph(Surv(dager, dod)~beh, data = pbc)
fit = survfit(Surv(dager, dod)~beh, data = pbc)
plot(cox.zph(cox))
plot(cox.zph(cox, transform = log))
plot(fit,fun="cloglog",lty=1:2,xlab="time",ylab="log(cum.haz)")
legend(50, 0, c("treatment","placebo"),lty=1:2,bty="n")
grid()

cox_res_uni = function(cox, cov_name){
  cs = summary(cox)
  pars = cs$conf.int
  sctest = cs$sctest[3]
  zph = cox.zph(cox)$table[3]
  zph_log = cox.zph(cox, transform = log)$table[3]
  
  #cov_name = 'treatment'
  df = data.frame(matrix(c(pars, sctest, zph, zph_log), nrow=1))
  names(df) = c(colnames(pars), 'score p', 'zph km p', 'zph log p')
  row.names(df) = cov_name
  df
}

library(xtable)
df = cox_res_uni(cox, 'Treatment')
caption = 'Cox regression on treatment covariate. score is the score test, zph is the proportional hazard test using scaled Schoenfeld-residuals.'
label = 'tab:cox_treatment'
tab =  xtable(df, caption = caption, label = label, 
              display = rep('g', 8), digits = 3)
print(tab, type = 'latex',  file = '~/stk4080/latex/tables/cox_treatment.tex')

#-----------
# Age
library(gam)
#cox = coxph(Surv(dager, dod)~factor(ald_groups), data = pbc)
cox = coxph(Surv(dager, dod)~ald, data = pbc)
summary(cox)
summary(cox)$sctest
martres = cox$residuals
plot(pbc$ald, martres)
ald = pbc$ald
printfig('cox_age_mart', height = 4)
plot(gam(martres~s(ald)), se=T, ylim=c(min(martres), max(martres)), 
     xlab='age', ylab='martingale residual')
points(ald, martres)
abline(0, 0, lty=4)
dev.off()

df = cox_res_uni(cox, 'Age')
caption = 'Cox regression on the age covariate. score is the score test, zph is the proportional hazard test using scaled Schoenfeld-residuals.'
label = 'tab:cox_age'
tab =  xtable(df, caption = caption, label = label, 
              display = rep('g', 8), digits = 3)
print(tab, type = 'latex',  file = '~/stk4080/latex/tables/cox_age.tex')

# Schoenfeld residual test for proportionality (lec 9)
# Small p-values indicate violation of proportionality assumption
cox.zph(cox)
plot(cox.zph(cox))
termplot(cox, se=TRUE) # can be used when we fit spline in cox regression

#---------
# Gender
cox = coxph(Surv(dager, dod)~factor(kjonn), data = pbc)
df = cox_res_uni(cox, 'Gender')
caption = 'Cox regression on the gender covariate. score is the score test, zph is the proportional hazard test using scaled Schoenfeld-residuals.'
label = 'tab:cox_gender'
tab =  xtable(df, caption = caption, label = label, 
              display = rep('g', 8), digits = 3)
print(tab, type = 'latex',  file = '~/stk4080/latex/tables/cox_gender.tex')

#-------
# Bilirubin
cox = coxph(Surv(dager, dod)~bil, data = pbc)
cox
martres = cox$residuals
bil = pbc$bil
printfig('cox_bil_mart', height = 4)
plot(gam(martres~s(bil)), se=T, ylim=c(min(martres), max(martres)), 
     xlab='bilirubin', ylab='martingale residual')
points(bil, martres)
abline(0, 0, lty=4)
dev.off()

df = cox_res_uni(cox, 'Bilirubin')
caption = 'Cox regression on the bilirubin covariate. score is the score test, zph is the proportional hazard test using scaled Schoenfeld-residuals.'
label = 'tab:cox_bilirubin'
tab =  xtable(df, caption = caption, label = label, 
              display = rep('g', 8), digits = 3)
print(tab, type = 'latex',  file = '~/stk4080/latex/tables/cox_bilirubin.tex')

cox = coxph(Surv(dager, dod)~pspline(bil), data = pbc)
cox
printfig('cox_bil_termplot', height = 4)
termplot(cox, se=TRUE, log='x', xlab='bilirubin', 
         ylab='partial for pspline(bilirubin)')
dev.off()
?termplot

cox = coxph(Surv(dager, dod)~log(bil), data = pbc)
cox
martres = cox$residuals
bil = pbc$bil
printfig('cox_bil_mart_log', height = 4)
plot(gam(martres~s(bil)), se=T, ylim=c(min(martres), max(martres)), 
     xlab='bilirubin', ylab='martingale residual', log='')
points(bil, martres)
abline(0, 0, lty=4)
dev.off()

df = cox_res_uni(cox, 'Bilirubin')
caption = 'Cox regression on the log-transformed bilirubin covariate. score is the score test, zph is the proportional hazard test using scaled Schoenfeld-residuals.'
label = 'tab:cox_bilirubin_log'
tab =  xtable(df, caption = caption, label = label, 
              display = rep('g', 8), digits = 3)
print(tab, type = 'latex',  file = '~/stk4080/latex/tables/cox_bilirubin_log.tex')

cox = coxph(Surv(dager, dod)~pspline(log(bil)), data = pbc)
cox


termplot(cox, se=TRUE)
plot(gam(martres~s(log(bil))), se=T, ylim=c(min(martres), max(martres)), 
     xlab='bilirubin', ylab='martingale residual')
points(bil, martres)
abline(0, 0, lty=4)
par(mfrow=c(1, 1))

#-------
# Albumin
cox = coxph(Surv(dager, dod)~alb, data = pbc)
cox
martres = cox$residuals
alb = pbc$alb
printfig('cox_alb_mart', height = 4)
plot(gam(martres~s(alb)), se=T, ylim=c(min(martres), max(martres)), 
     xlab='albumin', ylab='martingale residual')
points(alb, martres)
abline(0, 0, lty=4)
dev.off()

df = cox_res_uni(cox, 'Albumin')
caption = 'Cox regression on the albumin covariate. score is the score test, zph is the proportional hazard test using scaled Schoenfeld-residuals.'
label = 'tab:cox_albumin'
tab =  xtable(df, caption = caption, label = label, 
              display = rep('g', 8), digits = 3)
print(tab, type = 'latex',  file = '~/stk4080/latex/tables/cox_albumin.tex')

#############################################################
# c) Multivariat regresjon:
# Endelig utfører du en multivariat Cox-regresjonsanalyse hvor betydningen 
# av kovariatene studeres simultant. I denne analysen inng ̊ar det blant
# annet aa:
# (i) avgjøre hvilke kovariater som er av betydning for dødeligheten,
# (ii) avgjøre om effekten av behandling (eller andre kovariater) 
#      avhenger av verdien til øvrige kovariater (interaksjon),
# (iii) kontrollere modellens forutsetninger.

cox = coxph(Surv(dager, dod)~., data = pbc)
cox = coxph(Surv(dager, dod)~factor(beh) + ald + factor(kjonn) + log(bil) + alb, data = pbc)
cox
cox = coxph(Surv(dager, dod)~ald + factor(kjonn) + log(bil) + alb, data = pbc)
cox
cox = coxph(Surv(dager, dod)~ald + log(bil) + alb, data = pbc)
cox

library(lattice)
printfig('cor_heat', height = 5, width = 5)
x = cor(pbc[c('beh', 'ald', 'kjonn', 'bil', 'alb')])
row.names(x) = c('treat', 'age', 'gender', 'bil', 'alb')
colnames(x) = row.names(x)
levelplot(x, pretty=TRUE, xlab='', ylab='')
#levelplot(cor(pbc[c('beh', 'ald', 'kjonn', 'bil', 'alb')]), pretty=TRUE, xlab='',
#          ylab='')
dev.off()
?levelplot

plot(pbc$bil, pbc$alb, log = 'x', xlab = 'billirubin', ylab='albumin')
plot(pbc[c('beh', 'ald', 'kjonn', 'bil', 'alb')])

# Test first-order combinations without gender and treatment
cox = coxph(Surv(dager, dod)~ald * log(bil) + alb, data = pbc)
cox #not
cox = coxph(Surv(dager, dod)~ald * alb + log(bil), data = pbc)
cox #not
cox = coxph(Surv(dager, dod)~ald + log(bil) * alb, data = pbc)
cox #not

# Test all possible first-order interactions (without treatment)
cox = coxph(Surv(dager, dod)~ald*factor(kjonn) + log(bil) + alb, data = pbc)
cox # not
cox = coxph(Surv(dager, dod)~ald*log(bil) + factor(kjonn) + alb, data = pbc)
cox # ... not
cox = coxph(Surv(dager, dod)~ald*alb + factor(kjonn) + log(bil), data = pbc)
cox # not
cox = coxph(Surv(dager, dod)~ald + factor(kjonn)*log(bil) + alb, data = pbc)
cox 
cox = coxph(Surv(dager, dod)~ald + factor(kjonn)*alb + log(bil), data = pbc)
cox # not
cox = coxph(Surv(dager, dod)~ald + factor(kjonn) + log(bil)*alb, data = pbc)
cox # not
# Investigate the interesting without gender 
cox = coxph(Surv(dager, dod)~ald*log(bil) + alb, data = pbc)
cox # not
cox = coxph(Surv(dager, dod)~ald + (kjonn)*log(bil) + alb, data = pbc)
cox # probably because gender is so scewed. So it is basically just bil


# Test interactions with treatment
cox = coxph(Surv(dager, dod)~factor(beh)*ald + factor(kjonn) + log(bil) + alb, data = pbc)
cox 
cox = coxph(Surv(dager, dod)~factor(beh)*factor(kjonn) + ald + log(bil) + alb, data = pbc)
cox 
cox = coxph(Surv(dager, dod)~factor(beh)*log(bil) + ald + factor(kjonn) + alb, data = pbc)
cox 
cox = coxph(Surv(dager, dod)~factor(beh)*alb + ald + factor(kjonn) + log(bil), data = pbc)
cox 

cox = coxph(Surv(dager, dod)~factor(beh)*ald + log(bil) + alb, data = pbc)
cox 
cox = coxph(Surv(dager, dod)~factor(beh)*log(bil) + ald + alb, data = pbc)
cox 
cox = coxph(Surv(dager, dod)~factor(beh)*alb + ald + log(bil), data = pbc)
cox 
# ... none of them are significant


# Test model assumptions
cox = coxph(Surv(dager, dod)~ald + log(bil) + alb, data = pbc)
cox 
summary(cox)

4/10
printfig('cox_full_linearity', height = 3, width = 7.5)
par(mfrow=c(1, 3))
martres = cox$residuals
ald = pbc$ald
plot(gam(martres~s(ald)), se=T, ylim=c(min(martres), max(martres)), 
     xlab='age', ylab='martingale residual')
points(ald, martres)
abline(0, 0, lty=4)


martres = cox$residuals
bil = pbc$bil
plot(gam(martres~s(bil)), se=T, ylim=c(min(martres), max(martres)), 
     xlab='bil', ylab='martingale residual')
points(bil, martres)
abline(0, 0, lty=4)

martres = cox$residuals
alb = pbc$alb
plot(gam(martres~s(alb)), se=T, ylim=c(min(martres), max(martres)), 
     xlab='alb', ylab='martingale residual')
points(alb, martres)
abline(0, 0, lty=4)
dev.off()

cox.zph(cox)
cox.zph(cox, transform = log)

# Interpret model
cs = summary(cox)
cs
df = data.frame(cs$coefficients, check.names = FALSE)
row.names(df) = c('age', 'log(bilirubin)', 'albimin')
caption = 'Coefficinets for cox regression.'
label = 'tab:cox_final'
tab =  xtable(df, caption = caption, label = label, 
              display = rep('g', 6), digits = 3)
print(tab, type = 'latex',  file = '~/stk4080/latex/tables/cox_final.tex')

pars = cs$conf.int
pars
df = data.frame(pars, check.names = FALSE)
row.names(df) = c('age', 'log(bilirubin)', 'albimin')
caption = 'CI for cox regression.'
label = 'tab:cox_final_ci'
tab =  xtable(df, caption = caption, label = label, 
              display = rep('g', 5), digits = 3)
print(tab, type = 'latex',  file = '~/stk4080/latex/tables/cox_final_ci.tex')

cs$sctest
zph = cox.zph(cox)
df = data.frame(zph$table, check.names = FALSE)
row.names(df) = c('age', 'log(bilirubin)', 'albimin', 'GLOBAL')
caption = 'Test for proportional hazard using Schoenfeld-residuals. K-M transfromed.'
label = 'tab:cox_time_dep'
tab =  xtable(df, caption = caption, label = label, 
              display = rep('g', 4), digits = 3)
print(tab, type = 'latex',  file = '~/stk4080/latex/tables/cox_time_dep.tex')

zph_log = cox.zph(cox, transform = log)
df = data.frame(zph_log$table, check.names = FALSE)
row.names(df) = c('age', 'log(bilirubin)', 'albimin', 'GLOBAL')
caption = 'Test for proportional hazard using Schoenfeld-residuals. Log transformed.'
label = 'tab:cox_time_dep_log'
tab =  xtable(df, caption = caption, label = label, 
              display = rep('g', 4), digits = 3)
print(tab, type = 'latex',  file = '~/stk4080/latex/tables/cox_time_dep_log.tex')




cox = coxph(Surv(dager, dod)~pspline(ald) + factor(kjonn) + pspline(bil) + pspline(alb), data = pbc)
cox 
cox = coxph(Surv(dager, dod)~ald + factor(kjonn) + pspline(bil) + alb, data = pbc)
summary(cox)
cox = coxph(Surv(dager, dod)~ald + factor(kjonn) + log(bil) + alb, data = pbc)
summary(cox)
cox = coxph(Surv(dager, dod)~ald + log(bil) + alb, data = pbc)
summary(cox)
cox = coxph(Surv(dager, dod)~ald + factor(kjonn) + bil + alb, data = pbc)
summary(cox)
par(mfrow=c(2, 3), mai=c(1, 1, 0.1, 0.1))
termplot(cox, se=TRUE)
par(mfrow=c(1, 1))

