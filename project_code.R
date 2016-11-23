
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
 

# Not unique death / censor times
# We have tied data because of rounding to nearest day.
length(pbc$dager) 
length(unique(pbc$dager))
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
# b)
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
label = 'tab:cox_treatment'
tab =  xtable(df, caption = caption, label = label, 
              display = rep('g', 8), digits = 3)
print(tab, type = 'latex',  file = '~/stk4080/latex/tables/cox_age.tex')

# Schoenfeld residual test for proportionality (lec 9)
# Small p-values indicate violation of proportionality assumption
cox.zph(cox)
plot(cox.zph(cox))

?cox.zph
?termplot
termplot(cox, se=TRUE) # can be used when we fit spline in cox regression

cox = coxph(Surv(dager, dod)~kjonn, data = pbc)
?survfit
cox
cox = coxph(Surv(dager, dod)~bil_groups, data = pbc)
cox
cox = coxph(Surv(dager, dod)~alb_groups, data = pbc)
cox






