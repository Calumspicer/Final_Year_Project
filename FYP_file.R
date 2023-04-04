library("demography")
library("viridis")
library("tidyverse")
install.packages("bannerCommenter")
library(bannerCommenter)
#run the following command which saves the comment to clipboard
#then just paste on next line of script
banner("", emph = FALSE)  

# READ IN DATA
dat=read.demogdata(file="/Users/calum/Downloads/JPN/STATS/Mx_1x1.txt",popfile="/Users/calum/Downloads/JPN/STATS/Exposures_1x1.txt", type="mortality",label="Japan",skip=0,scale=1 )
dat2=read.demogdata(file="/Users/calum/Downloads/Netherlands_Mx_1x1.txt",popfile="/Users/calum/Downloads/Exposures_1x1.txt", type="mortality",label="Netherlands",skip=0,scale=1 )
dat_sub=extract.years(dat,c(1950:2000))
dat_sub_NED=extract.years(dat2,c(1950:2000))
dat_pred=extract.years(dat,c(2001:2020))
plot(dat_sub,series="male")
dat_sub65=extract.ages(dat_sub,65,TRUE)

dat_sub$rate[1]

#PERFORM LCA
#plot(dat_sub)
lca_JAP=lca(dat_sub,series=names(dat_sub$rate)[1],years=dat_sub$year,ages=dat_sub$age,max.age=95,interpolate=TRUE)
#plot(lca_JAP)
lca_NED=lca(dat_sub_NED,series=names(dat_sub_NED$rate)[1],years=dat_sub_NED$year,ages=dat_sub_NED$age,max.age=95,interpolate=TRUE)


#CREATE DATA FRAMES OF EACH PARAMETER
japax=data.frame(lca_JAP$ax,lca_JAP$age)
nedax=data.frame(lca_NED$ax,lca_JAP$age)
japbx=data.frame(lca_JAP$bx,lca_JAP$age)
nedbx=data.frame(lca_NED$bx,lca_JAP$age)
japkt=data.frame(lca_JAP$kt,lca_JAP$year)
nedkt=data.frame(lca_NED$kt,lca_JAP$year)


#plot(lca_NED)

#sum(lca_JAP$kt)



########################

#ARIMA MODELLING
auto.arima(lca_JAP$kt)
sqrt(12.23)

fit=Arima(lca_JAP$kt,order=c(0,1,1),include.drift = TRUE)

predic=forecast(fit,h=30)
kpred=predic$mean
##########################

#CALCULATE DIFFERENCES
diff=1:10
for(i in 2:length(lca_JAP$kt)){
  diff[i]=lca_JAP$kt[i]-lca_JAP$kt[i-1]
}
kt_sd=sd(diff[-1])





#japkt2_sim=data.frame(k_x,lca_JAP$year)

k_x=1:10
m_x=1:10
k_x[1]=lca_JAP$kt[1]
#typeof(k_x[1])

theta=(tail(lca_JAP$kt,1)-k_x[1])/(length(lca_JAP$kt)-1)



for(i in 2:51){
  k_x[i]=k_x[i-1]+theta+rnorm(1,mean=0,sd=kt_sd)
}

k_sample=function(){
  for(i in 2:51){
    k_x[i]=k_x[i-1]+theta+rnorm(1,mean=0,sd=kt_sd)
  }
  return(k_x)
}


N=10
M=51
blah <- data.frame(matrix(c(1950:2000), nrow=M, ncol=N))
for (i in 2:N) {
  yourResults <- k_sample()
  blah[,i] <- yourResults
  yt2=c(yt2,yourResults)
}





mort=1:30
kpred=1:30
kpred[1]=tail(japkt$lca_JAP.kt,n=1)


for(i in 1:30){
  kpred[i+1]=kpred[i]+theta+rnorm(1,mean=0,sd=kt_sd)
  mort[i]=japax$lca_JAP.ax[65+i]+  japbx$lca_JAP.bx[65+i]*kpred[i]
}

japax$lca_JAP.ax[66]
japbx$lca_JAP.bx[66]
kpred[1]

sum(japbx$lca_JAP.bx)
sum(japkt$lca_JAP.kt)

mort
length(2000:2030)
dat$rate[1]$female[66,]

px=1-exp(-exp(mort))
#matrix(data=c(2000:2029,exp(mort),exp(-exp(mort)),1-exp(-exp(mort))),nrow=4)
cohort=matrix(data=c(2001:2030,66:95,exp(mort),exp(-exp(mort)),1-exp(-exp(mort))),ncol=5)
colnames(cohort)=c("Year","Age","Mx,t","px","qx")
cohort

px[1]
#WANG ADJUSTED
wang=pnorm(qnorm(px,mean=0,sd=1)-0.05,mean=0,sd=1)
wang2=pnorm(qnorm(px,mean=0,sd=1)-0.1,mean=0,sd=1)
pxts=ts(px,start=2001)
wangts=ts(wang,start=2001)
wang2ts=ts(wang2,start=2001)
plot(pxts,ylab="qx values")
lines(wangts,col="red")
lines(wang2ts,col="blue")
legend(2000,0.1, legend=c("qx", "0.05 (Wang)","0.1 (Wang)"),
       col=c("black", "red","blue"), lty=1, cex=0.8)



#OUTPUT NEW COHORT TABLE
px_adj=1-wang
mx_adj=log(1/px_adj)

remain=1:10
remain[1]=1
for(i in 2:31){
remain[i]=remain[i-1]*px_adj[i-1]  
}

remain

cohort_adj=matrix(data=c(2001:2030,66:95,mx_adj,px_adj,wang,remain[-1]),ncol=6)
colnames(cohort_adj)=c("Year","Age","Mx,t","px","qx","Remaining")
cohort_adj

px_adj
log(1/px_adj)
mort
mx_adj




exp(-mx_adj)





#https://www.ecb.europa.eu/stats/financial_markets_and_interest_rates/euro_area_yield_curves/html/index.en.html

intRate=0.02159

discFact=1/(1+intRate)
#######################
#FIXED LEG

#SUM(SURVIVAL_PROB*DISCOUNT_FACTOR)
fix=1:30
for(i in 1:30){
fix[i]=px_adj[i]*(discFact^i)  
}


fixedLeg=sum(fix)



#FLOATING LEG

#SUM(SURVIVAL_PROB*1_YEAR_COHORT_MORTALITY_RATE*DISCOUNT_FACTOR)
float=1:30
for(i in 1:30){
  float[i]=px_adj[i]*mx_adj[i]*(discFact^i)  
}
float

floatingLeg=sum(float)


########################

#Longevity Swap Rate

# floating_leg/fixed_leg

longRate=floatingLeg/fixedLeg


cc=forecast(lca_JAP)

cc$rate$female[,1]















##################################################################
##################################################################



##################################################################
##################################################################
#############################################
#CLEANED UP



lca_JAPf=lca(dat_sub,series=names(dat_sub$rate)[1],years=dat_sub$year,ages=dat_sub$age,max.age=95,interpolate=FALSE)
lca_JAPm=lca(dat_sub,series=names(dat_sub$rate)[2],years=dat_sub$year,ages=dat_sub$age,max.age=95,interpolate=FALSE)

lca_NEDf=lca(dat_sub_NED,series=names(dat_sub_NED$rate)[1],years=dat_sub_NED$year,ages=dat_sub_NED$age,max.age=95,interpolate=FALSE)
lca_NEDm=lca(dat_sub_NED,series=names(dat_sub_NED$rate)[2],years=dat_sub_NED$year,ages=dat_sub_NED$age,max.age=95,interpolate=FALSE)

diffJapF=1:10
diffJapM=1:10
diffNedF=1:10
diffNedM=1:10


for(i in 2:length(lca_JAP$kt)){
  diffJapF[i-1]=lca_JAPf$kt[i]-lca_JAPf$kt[i-1]
  diffJapM[i-1]=lca_JAPm$kt[i]-lca_JAPm$kt[i-1]
  diffNedF[i-1]=lca_NEDf$kt[i]-lca_NEDf$kt[i-1]
  diffNedM[i-1]=lca_NEDm$kt[i]-lca_NEDm$kt[i-1]
}

ktsd_JapF=sd(diffJapF)
ktsd_JapM=sd(diffJapM)
ktsd_NedF=sd(diffNedF)
ktsd_NedM=sd(diffNedM)


theta_JapF=mean(diffJapF)
theta_JapM=mean(diffJapM)
theta_NedF=mean(diffNedF)
theta_NedM=mean(diffNedM)


ktJapF=1:10
ktJapM=1:10
ktNedF=1:10
ktNedM=1:10


mortJapF=1:10
mortJapM=1:10
mortNedF=1:10
mortNedM=1:10


for(i in 1:50){
  ktJapF[i]=tail(lca_JAPf$kt,1)+theta_JapF*i
  ktJapM[i]=tail(lca_JAPm$kt,1)+theta_JapM*i
  ktNedF[i]=tail(lca_NEDf$kt,1)+theta_NedF*i
  ktNedM[i]=tail(lca_NEDm$kt,1)+theta_NedM*i
}
  
for(i in 1:30){  
  mortJapF[i]=lca_JAPf$ax[65+i]+lca_JAPf$bx[65+i]*ktJapF[i]
  mortJapM[i]=lca_JAPm$ax[65+i]+lca_JAPm$bx[65+i]*ktJapM[i]
  mortNedF[i]=lca_NEDf$ax[65+i]+lca_NEDf$bx[65+i]*ktNedF[i]
  mortNedM[i]=lca_NEDm$ax[65+i]+lca_NEDm$bx[65+i]*ktNedM[i]
}


compare.demogdata(dat_pred,forecast(lca_JAPf,jumpchoice = "actual"))

mxJapF=exp(mortJapF)
mxJapM=exp(mortJapM)
mxNedF=exp(mortNedF)
mxNedM=exp(mortNedM)

pxJapF=exp(-exp(mortJapF))
pxJapM=exp(-exp(mortJapM))
pxNedF=exp(-exp(mortNedF))
pxNedM=exp(-exp(mortNedM))


qxJapF=1-exp(-exp(mortJapF))
qxJapM=1-exp(-exp(mortJapM))
qxNedF=1-exp(-exp(mortNedF))
qxNedM=1-exp(-exp(mortNedM))

##################################################################
##################################################################
#WANG TRANSFORM

wangTrans=function(w){
  a=pnorm(qnorm(1-exp(-exp(mortJapF)),mean=0,sd=1)-w,mean=0,sd=1)
  b=pnorm(qnorm(1-exp(-exp(mortJapM)),mean=0,sd=1)-w,mean=0,sd=1)
  c=pnorm(qnorm(1-exp(-exp(mortNedF)),mean=0,sd=1)-w,mean=0,sd=1)
  d=pnorm(qnorm(1-exp(-exp(mortNedM)),mean=0,sd=1)-w,mean=0,sd=1)
  
  ab=matrix(data=c(a,b,c,d),ncol=4)
  colnames(ab)=c("JapF","JapM","NedF","NedM")
  return(ab)
}

wangTrans(0.01)


pxadj=1-wangTrans(0.01)


mxadj=log(1/pxadj)

log(mxadj)

lf=1:10
lf2=1:10

life_expect_func=function(c,g,w){
  lf=1:10
  lf2=1:10
  pxa=1-wangTrans(w)
for(i in 1:30){
lf[i]=sum(tail(pxa[,(c+g-1)],31-i))  
lf2[i]=lf[i]+65+i
}
return(ts(lf2,start=2001))
}
lfts=ts(lf2,start=2001)

plot(lfts)





##################################################################
##################################################################

#COHORT TABLE

cohort_mat=function(c,g,w){
  a=wangTrans(w)
b=a[,c+g-1]
c=exp(-b)
d=1-c
t=1:3
for(i in 1:29){
  t[i+1]=t[i]*c[i]
}
ab=matrix(data=c(2001:2030,66:95,b,c,d,t),ncol=6)
colnames(ab)=c("Year","Age","Exp.Mort","px","qx","Remaining")
return(ab)
}


cohort_mat(1,1,0)

intRate=0.02159

discFact=1/(1+intRate)
#######################
#FIXED LEG

#SUM(SURVIVAL_PROB*DISCOUNT_FACTOR)
fixedLeg_calc=function(c,g,w,r){
fix=1:30
pp=cohort_mat(c,g,w)[,4]
discFact=1/(1+r)
for(i in 1:30){
  fix[i]=pp[i]*(discFact^i)  
}


fixedLeg=sum(fix)
return(fixedLeg)
}

fixedLeg_calc(1,1,0.1,intRate)

#FLOATING LEG

#SUM(SURVIVAL_PROB*1_YEAR_COHORT_MORTALITY_RATE*DISCOUNT_FACTOR)
float=1:30
for(i in 1:30){
  float[i]=px_adj[i]*mx_adj[i]*(discFact^i)  
}
float

floatingLeg=sum(float)

floatingLeg_calc=function(c,g,w,r){
  float=1:30
  pp=cohort_mat(c,g,w)[,4]
  mm=cohort_mat(c,g,w)[,3]
  discFact=1/(1+r)
  for(i in 1:30){
    float[i]=pp[i]*mm[i]*(discFact^i)  
  }
  
  
  floatingLeg=sum(float)
  return(floatingLeg)
}

floatingLeg_calc(1,1,0.1,intRate)
########################

#Longevity Swap Rate

# floating_leg/fixed_leg

longRate=floatingLeg/fixedLeg

###############################

#PLOTs

###############################

lcaJap=function(g){
  lca_JAP2=lca(dat_sub,series=names(dat_sub$rate)[g],years=dat_sub$year,ages=dat_sub$age,max.age=95,interpolate=FALSE)
return(lca_JAP2)
  }

#LCA
kt_sum=1:4
for(i in 1:2){
lca_JAP=lca(dat_sub,series=names(dat_sub$rate)[i],years=dat_sub$year,ages=dat_sub$age,max.age=95,interpolate=FALSE)
kt_sum[i]=sum(lca_JAP$kt)
plot(lca_JAP)
lca_NED=lca(dat_sub_NED,series=names(dat_sub_NED$rate)[i],years=dat_sub_NED$year,ages=dat_sub_NED$age,max.age=95,interpolate=FALSE)
plot(lca_NED)
kt_sum[2+i]=sum(lca_NED$kt)
}
kt_sum



##################################################################

#DEATH RATES PLOT

par(mfrow=c(1,2))

#JAPAN PLOTS
plot(dat_sub,series="female",main="Female Death Rates (1950-2000)")
plot(dat_sub,series="male",main="Male Death Rates (1950-2000)")

#NETHERLANDS PLOTS
plot(dat,series="female",main="Female Death Rates (1950-2000)")
plot(dat,series="male",main="Male Death Rates (1950-2000)")



##################################################################

#LCA PLOTS
par(mfrow=c(1,2))





plot(ts(lca_JAPf$ax,start=0),main=expression("A"[x]*"Values for Females"),ylab=expression("A"[x]*"Values"),xlab="Age")
lines(ts(lca_NEDf$ax,start=0),col="green")
legend(0,-2, legend=c("Japan","Netherlands"),
       col=c("black","green"), lty=1, cex=0.8)


plot(ts(lca_JAPm$ax,start=0),main=expression("A"[x]*"Values for Males"),ylab=expression("A"[x]*"Values"),xlab="Age")
lines(ts(lca_NEDm$ax,start=0),col="green")


plot(ts(lca_JAPf$bx,start=0),main=expression("B"[x]*"Values for Females"),ylab=expression("B"[x]*"Values"),xlab="Age",ylim=c(0.002,0.04))
lines(ts(lca_NEDf$bx,start=0),col="green")
legend(40,0.038, legend=c("Japan","Netherlands"),
       col=c("black","green"), lty=1, cex=0.8)

plot(ts(lca_JAPm$bx,start=0),main=expression("B"[x]*"Values for Males"),ylab=expression("B"[x]*"Values"),xlab="Age",ylim=c(0.002,0.04))
lines(ts(lca_NEDm$bx,start=0),col="green")

plot(ts(lca_JAPf$kt,start=1950),main=expression("K"[t]*"Values for Females"),ylab=expression("K"[t]*"Values"),xlab="Year",ylim=c(-100,100))
lines(ts(lca_NEDf$kt,start=1950),col="green")
legend(1970,90, legend=c("Japan","Netherlands"),
       col=c("black","green"), lty=1, cex=0.8)



plot(ts(lca_JAPm$kt,start=1950),main=expression("K"[t]*"Values for Males"),ylab=expression("K"[t]*"Values"),xlab="Year",ylim=c(-100,100))
lines(ts(lca_NEDm$kt,start=1950),col="green")





exp(lca_JAPf$ax[21])
exp(lca_NEDf$ax[21])




plot(lca_JAPm)
plot(lca_NEDf)
plot(lca_NEDm)


##################################################################
##################################################################

#WANG PLOTS

par(mfrow=c(1,2))

wangts=ts(wangTrans(0.0),start=2001)
wangts2=ts(wangTrans(0.005),start=2001)
wangts3=ts(wangTrans(0.05),start=2001)
wangts4=ts(wangTrans(0.1),start=2001)
wangts5=ts(wangTrans(0.2),start=2001)
wangts6=ts(wangTrans(0.3),start=2001)

wang2ts=ts(wang2,start=2001)
plot(wangts[,4],ylab="qx values",ylim=c(0,0.2),main="Netherlands Males")
lines(wangts3[,4],col="red")
lines(wangts4[,4],col="blue")
lines(wangts5[,4],col="green")
lines(wangts6[,4],col="purple")
legend(2000,0.15, legend=c("qx", "0.05 (Wang)","0.1 (Wang)","0.2 (Wang)","0.3 (Wang)"),
       col=c("black", "red","blue","green","purple"), lty=1, cex=0.8)



##################################################################
##################################################################


cohort_mat(1,1,0.1)


write_csv(data.frame(cohort_mat(1,1,0.1)),file="/Users/calum/Library/CloudStorage/OneDrive-Personal/FYP/cohort_table.csv")


##################################################################
##################################################################

#Life Expectancy
gen=2
coun=2
plot(life_expect_func(coun,gen,0),ylim=c(94,96),main="Male Life Expectancy",ylab="Life Expectancy",xlab="Year")
lines(life_expect_func(coun,gen,0.05),col="red")
lines(life_expect_func(coun,gen,0.1),col="blue")
lines(life_expect_func(coun,gen,0.2),col="green")
lines(life_expect_func(coun,gen,0.3),col="purple")
legend(2000,96, legend=c(expression(lambda*"=0"), expression(lambda*"=0.05"),expression(lambda*"=0.1"),expression(lambda*"=0.2"),expression(lambda*"=0.3")),
       col=c("black", "red","blue","green","purple"), lty=1, cex=0.8)



lifetable(dat_sub,series=names(dat_sub$rate)[1],years = 2000,max.age = 95)


lca

par(mfrow=c(1,1))
plot(ts(lca_JAPf$kt,start=1950),xlim=c(1950,2030),ylim=c(-200,100))
lines(ts(ktJapF,start=2001))

plot(ts(lca_JAPm$kt,start=1950),xlim=c(1950,2030),ylim=c(-200,100))
lines(ts(ktJapM,start=2001))

lca_JAPf$female[66,]

plot(ts(lca_JAPf$female[66,],start=1950),xlim=c(1950,2030),ylim=c(0,0.015))
lines(ts(1-exp(-exp(mortJJ)),start=2001))

print(c(ktsd_JapF,ktsd_JapM,ktsd_NedF,ktsd_NedM))

mortJJ=1:10
for(i in 1:30){
 
  
  mortJJ[i]=lca_JAPf$ax[66]+lca_JAPf$bx[66]*ktJapF[i]
 
}


plot(lca_JAPf$kt*lca_JAPf$bx[66],xlim=c(1950,2030),ylim=c(-2,1))
lines(ts(ab,start=2001))
ab=ktJapF*lca_JAPf$bx[66]
ab=1-exp(-exp(mortJJ))

plot(ts(lca_JAPf$female[66,],start=1950),xlim=c(1950,2030),ylim=c(0,0.015))
lines(ts(ab,start=2001))


1-exp(-exp(lca_JAPf$ax[66]+lca_JAPf$kt*lca_JAPf$bx[66]))
plot(ts(lca_JAPf$female[66,],start=1950),ylim=c(0,0.04))
lines(ts(1-exp(-exp(lca_JAPf$ax[66]+lca_JAPf$kt*lca_JAPf$bx[66])),start=1950))

plot(forecast(lca_JAPf))

lx=forecast(lca_JAPf,h=20)
lx2=lx$rate$female

plot(ts(lx2[66,],start=2001),ylim=c(0,0.006))

lines(ts(dat_pred$rate$female[66,],start=2001))
lines(ts(ab,start=2001))



test=lifetable(dat_sub,series=names(dat$rate)[1],years = 2000)
test

test$qx[66]

mortJJ


sum(tail(test$lx,n=35))/test$lx[66]
test$lx[66]

life_expect_func2=function(c,g,w){
  lf=1:10
  lf2=1:10
  pxa=1-wangTrans(w)
  for(i in 1:30){
    lf[i]=sum(tail(pxa[,(c+g-1)],31-i))  
    lf2[i]=lf[i]+65+i
  }
  return(ts(lf2,start=2001))
}

1-wangTrans(0)

lf_JAPf=1:5
lf_JAPm=1:5
lf_NEDf=1:5
lf_NEDm=1:5


lf_JAPf[1]=sum(1-wangTrans(0)[,1])
lf_JAPm[1]=sum(1-wangTrans(0)[,2])
lf_NEDf[1]=sum(1-wangTrans(0)[,3])
lf_NEDm[1]=sum(1-wangTrans(0)[,4])


life.expectancy(dat_sub,series=names(dat_sub$rate)[1],years=2000,age=65)


xtr=forecast(lca_JAPf,30)
xtr$rate$female[65,]

xtr$kt.f
psi=1:10
for(i in 1:29){
psi[i]=xtr$rate$female[65+i,i]  
}

sum(1-exp(-exp(psi)))

dat_pred$rate$female[66,]





dat_sub$rate$female[,50]
lca_JAPf$female[,50]
forecast()
################################################
#LIFE EXPECTANCY PROJECTIONS
px_r=1:30
for(i in 1:30){
  px_r[i]=test$Lx[65+i]
}
sum(px_r)
qx_r2000=1:30
px_r2000=1:30
lx_r2000=1:30


for(i in 1:30){
  qx_r2000[i]=test$qx[65+i]
  px_r2000[i]=1-qx_r[i]
  lx_r2000[i+1]=lx_r[i]*px_r[i]
}


for(i in 1:50){
  ktJapF[i]=tail(lca_JAPf$kt,1)+theta_JapF*i
  ktJapM[i]=tail(lca_JAPm$kt,1)+theta_JapM*i
  ktNedF[i]=tail(lca_NEDf$kt,1)+theta_NedF*i
  ktNedM[i]=tail(lca_NEDm$kt,1)+theta_NedM*i
}

for(i in 1:30){  
  mortJapF[i]=lca_JAPf$ax[65+i]+lca_JAPf$bx[65+i]*ktJapF[5]
  mortJapM[i]=lca_JAPm$ax[65+i]+lca_JAPm$bx[65+i]*ktJapM[5]
  mortNedF[i]=lca_NEDf$ax[65+i]+lca_NEDf$bx[65+i]*ktNedF[5]
  mortNedM[i]=lca_NEDm$ax[65+i]+lca_NEDm$bx[65+i]*ktNedM[5]
}



qx_r2005=1:30
px_r2005=1:30
lx_r2005=1:30
for(i in 1:30){
  px_r2005[i]=exp(-exp(mortJapF[i]))
  qx_r2005[i]=1-px_r2005[i]
  
  lx_r2005[i+1]=lx_r2005[i]*px_r2005[i]
}

sum(lx_r2005)



for(i in 1:30){  
  mortJapF[i]=lca_JAPf$ax[65+i]+lca_JAPf$bx[65+i]*ktJapF[10]
  mortJapM[i]=lca_JAPm$ax[65+i]+lca_JAPm$bx[65+i]*ktJapM[10]
  mortNedF[i]=lca_NEDf$ax[65+i]+lca_NEDf$bx[65+i]*ktNedF[10]
  mortNedM[i]=lca_NEDm$ax[65+i]+lca_NEDm$bx[65+i]*ktNedM[10]
}



qx_r2010=1:30
px_r2010=1:30
lx_r2010=1:30
for(i in 1:30){
  px_r2010[i]=exp(-exp(mortJapF[i]))
  qx_r2010[i]=1-px_r2010[i]
  
  lx_r2010[i+1]=lx_r2010[i]*px_r2010[i]
}

sum(lx_r2010)

abc=c()
lx_JapF=1:4
lx_JapM=1:4
lx_NedF=1:4
lx_NedM=1:4

for(j in 1:4){
  qx_rJapF=1:30
  qx_rJapM=1:30
  qx_rNedF=1:30
  qx_rNedM=1:30
  
  px_rJapF=1:30
  px_rJapM=1:30
  px_rNedF=1:30
  px_rNedM=1:30
  
  lx_rJapF=1:30
  lx_rJapM=1:30
  lx_rNedF=1:30
  lx_rNedM=1:30
  

  
  for(i in 1:30){
    mortJapF[i]=lca_JAPf$ax[65+i]+lca_JAPf$bx[65+i]*ktJapF[i+5*j]
    mortJapM[i]=lca_JAPm$ax[65+i]+lca_JAPm$bx[65+i]*ktJapM[i+5*j]
    mortNedF[i]=lca_NEDf$ax[65+i]+lca_NEDf$bx[65+i]*ktNedF[i+5*j]
    mortNedM[i]=lca_NEDm$ax[65+i]+lca_NEDm$bx[65+i]*ktNedM[i+5*j]
    px_rJapF[i]=exp(-exp(mortJapF[i]))
    qx_rJapF[i]=1-px_rJapF[i]
    
    lx_rJapF[i+1]=lx_rJapF[i]*px_rJapF[i]
    
    
    px_rJapM[i]=exp(-exp(mortJapM[i]))
    qx_rJapM[i]=1-px_rJapM[i]
    lx_rJapM[i+1]=lx_rJapM[i]*px_rJapM[i]
    
    
    px_rNedF[i]=exp(-exp(mortNedF[i]))
    qx_rNedF[i]=1-px_rNedF[i]
    lx_rNedF[i+1]=lx_rNedF[i]*px_rNedF[i]
    
    
    px_rNedM[i]=exp(-exp(mortNedM[i]))
    qx_rNedM[i]=1-px_rNedM[i]
    lx_rNedM[i+1]=lx_rNedM[i]*px_rNedM[i]
  }
  
  lx_JapF[j]=sum(lx_rJapF)
  lx_JapM[j]=sum(lx_rJapM)
  lx_NedF[j]=sum(lx_rNedF)
  lx_NedM[j]=sum(lx_rNedM)
  abc=c(abc,lx_rJapF)
  
}

abcd

lx_JapF
plot(ts(lx_JapF,start=2000,deltat=5),ylim=c(15,30),ylab=expression("e"[0]),main="Remaining Lifetime of a 65 y/o")
lines(ts(lx_JapM,start=2000,deltat=5),col="blue")
lines(ts(lx_NedF,start=2000,deltat=5),col="red")
lines(ts(lx_NedM,start=2000,deltat=5),col="green")
legend(2000,30, legend=c("Japan (f)", "Japan (m)", "Netherlands (f)", "Netherlands (m)"),
       col=c("black", "blue","red","green"), lty=1, cex=0.8)


lx_JapF[10]
abc
abcd=matrix(data=abc,ncol=10)
abcd











for(i in 1:30){
  mortJapF[i]=lca_JAPf$ax[65+i]+lca_JAPf$bx[65+i]*ktJapF[i]
  mortJapM[i]=lca_JAPm$ax[65+i]+lca_JAPm$bx[65+i]*ktJapM[i]
  mortNedF[i]=lca_NEDf$ax[65+i]+lca_NEDf$bx[65+i]*ktNedF[i]
  mortNedM[i]=lca_NEDm$ax[65+i]+lca_NEDm$bx[65+i]*ktNedM[i]
  px_rJapF[i]=exp(-exp(mortJapF[i]))
  qx_rJapF[i]=1-px_rJapF[i]
  
  lx_rJapF[i+1]=lx_rJapF[i]*px_rJapF[i]
  
  
  px_rJapM[i]=exp(-exp(mortJapM[i]))
  qx_rJapM[i]=1-px_rJapM[i]
  lx_rJapM[i+1]=lx_rJapM[i]*px_rJapM[i]
  
  
  px_rNedF[i]=exp(-exp(mortNedF[i]))
  qx_rNedF[i]=1-px_rNedF[i]
  lx_rNedF[i+1]=lx_rNedF[i]*px_rNedF[i]
  
  
  px_rNedM[i]=exp(-exp(mortNedM[i]))
  qx_rNedM[i]=1-px_rNedM[i]
  lx_rNedM[i+1]=lx_rNedM[i]*px_rNedM[i]
}

#https://www.ecb.europa.eu/stats/financial_markets_and_interest_rates/euro_area_yield_curves/html/index.en.html

intRate=0.02159

discFact=1/(1+intRate)
#######################
#FIXED LEG

#SUM(SURVIVAL_PROB*DISCOUNT_FACTOR)
fix=1:30

xp=30
ind=4
fix_mat=matrix(nrow=xp,ncol=4)
for(i in 1:xp){
  fix[i]=px_adj[i]*(discFact^i)
  fix_mat[i,1]=px_rJapF[i]*(discFact^i)
  fix_mat[i,2]=px_rJapM[i]*(discFact^i)
  fix_mat[i,3]=px_rNedF[i]*(discFact^i)
  fix_mat[i,4]=px_rNedM[i]*(discFact^i)
  
}

fixedLeg2=matrix(nrow=4,ncol=4)
fixedLeg2[,ind]=colSums(fix_mat)
fixedLeg



#FLOATING LEG
pxa=matrix(px_)
#SUM(SURVIVAL_PROB*1_YEAR_COHORT_MORTALITY_RATE*DISCOUNT_FACTOR)
float=1:30
float_mat=matrix(nrow=xp,ncol=4)
for(i in 1:xp){
  float[i]=px_adj[i]*mx_adj[i]*(discFact^i)  
  float_mat[i,1]=px_rJapF[i]*exp(mortJapF[i])*(discFact^i)
  float_mat[i,2]=px_rJapM[i]*exp(mortJapM[i])*(discFact^i)
  float_mat[i,3]=px_rNedF[i]*exp(mortNedF[i])*(discFact^i)
  float_mat[i,4]=px_rNedM[i]*exp(mortNedM[i])*(discFact^i)
}
float

float_mat
floatingLeg=colSums(float_mat)
floatingLeg2=matrix(nrow=4,ncol=4)
floatingLeg2[,ind]=colSums(float_mat)
floatingLeg

########################

#Longevity Swap Rate

# floating_leg/fixed_leg
longRate2=matrix(nrow=4,ncol=4)
longRate2[,ind]=colSums(float_mat)/colSums(fix_mat)
longRate
floatingLeg
fixedLeg


plot(c(5,10,20,30),fixedLeg2[1,],ylim=c(4,22),type="o",xlab="Contract Duration (Years)",ylab="Fixed Leg NPV",main="Fixed Leg Value for Different Contract Durations")
lines(c(5,10,20,30),fixedLeg2[2,],col="red",type="o",lwd=2)
lines(c(5,10,20,30),fixedLeg2[3,],col="blue",type="o")
lines(c(5,10,20,30),fixedLeg2[4,],col="green",type="o")
legend(5,22, legend=c("Japan (f)", "Japan (m)", "Netherlands (f)", "Netherlands (m)"),
       col=c("black" ,"red","blue","green"), lty=1, cex=0.8)

