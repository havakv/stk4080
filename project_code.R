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
# 
# Not unique death / censor times
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
fit = survfit(Surv(dager, dod)~beh, data = pbc)
plot(fit, lty=1:2, ylab='survival', xlab='days')
legend(1,0.2,c("treatment","placebo"), lty=1:2, bty='n')
grid()

library(dplyr)
df = tbl_df(pbc)
glimpse(df)
df %>% group_by(kjonn) %>% summarise(avg = mean(dager))

#----
nb_groups = 3
ald_groups = hist(pbc$ald, breaks = nb_groups)$breaks
pbc$ald_groups = 0
ald_groups[1:(length(ald_groups)-1)]
for (i in ald_groups[1:(length(ald_groups)-1)]) {
  pbc$ald_groups[pbc$ald > i] = i
}
gr = sort(unique(pbc$ald_groups))

fit = survfit(Surv(dager, dod)~ald_groups, data = pbc)
plot(fit, lty=1:length(gr), ylab='survival', xlab='days')
legend(1,0.2, gr, lty=1:length(gr), bty='n')
grid()

#----
fit = survfit(Surv(dager, dod)~kjonn, data = pbc)
plot(fit, lty=1:2, ylab='survival', xlab='days')
legend(1,0.2,c("male","female"), lty=1:2, bty='n')
grid()

#----
pbc$bil_groups = 0
hist(pbc$bil, breaks = 20)
gr = c(0, 1, 3, 7, 15)
for (i in gr) {
  pbc$bil_groups[pbc$bil > i] = i
}

df %>% group_by(as.integer(pbc$bil_groups)) %>% summarise(length(dager))

fit = survfit(Surv(dager, dod)~bil_groups, data = pbc)
plot(fit, lty=1:length(gr), ylab='survival', xlab='days')
legend(1,0.5, gr, lty=1:length(gr), bty='n')
grid()

#----
pbc$alb_groups = 0
hist(pbc$alb, breaks = 20)
gr = c(0, 3, 3.5, 4)
for (i in gr) {
  pbc$alb_groups[pbc$alb > i] = i
}

fit = survfit(Surv(dager, dod)~alb_groups, data = pbc)
plot(fit, lty=1:length(gr), ylab='survival', xlab='days')
legend(1, 0.4, gr, lty=1:length(gr), bty='n')
grid()
