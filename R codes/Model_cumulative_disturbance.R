###### Bayesian model to quantify the influence of enviromental impacts on coral cover measurements
# from in Vercelloni et al. Disentangling the impacts of disturbance on coral reefs, a modern approach to an old problem 

rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(rjags)
library(runjags)
library(tidymcmc) # https://github.com/bragqut/tidymcmc


# Read the data from downloaded github folder  
df<-load("../Data/Cumulative_disturbance.Rdata")

##################################################################################### Data manipulation 
##Keep corals above Cairns

df_coral<-df%>%
  filter(lat>-17& lng<147)%>%
  arrange(stid,year)

## avoid null value

df_coral<-df_coral%>%
  mutate(cover= ifelse(cover==0,0.001,cover))%>%
  filter(!is.na(cover))

## Keep sub-transect repeated 3 times 
df_coral_full<-df_coral%>%group_by(stid)%>%tally()%>%filter(n>2)
df_coral<-inner_join(df_coral,df_coral_full)

##################################################################################### The model 
model<-  "
model{

# likelihood for beta distribution
for(i in 1:n) { 
y[i] ~ dbeta(alpha[i], beta[i]) 
alpha[i] <- mu[i] * phi
beta[i]  <- (1-mu[i]) * phi
logit(mu[i]) <- b4[Group[i],StID[i]]*Ita.Bleaching[i]+b5[Group[i],StID[i]]*Ita.Nathan.Bleaching[i]

#Assess model fit using Pearson residuals


alta[i] <- alpha[i] + beta[i]
expected[i] <- alpha[i] / alta[i] 
variance[i] <- alpha[i] * beta[i] / alta[i]*alta[i]*(alta[i]+1)
residual[i] <- (y[i]-expected[i]) / sqrt(variance[i])                                     
sq[i] <- residual[i]^2 

}

# priors

for (h in 1:ngroup){
for (j in 1:nstid){
b4[h,j]~dnorm(b.hat4[h], taub4[j])
b5[h,j]~dnorm(b.hat5[h], taub5[j])
}}

for (h in 1:ngroup){
b.hat4[h]~dnorm(b.hat4.0,taub.0)
b.hat5[h]~dnorm(b.hat5.0,taub.0)
}


phi ~ dgamma(.1,.1)
b.hat4.0~dnorm(0,1.0E-2)
b.hat5.0~dnorm(0,1.0E-2)
taub.0~dgamma(1.0E-2, 1.0E-2)
}"


##Inits

inits1 <- list(b.hat1=rnorm(1, mean = 0, sd = 1.0E-2), b.hat2=rnorm(1, mean = 0, sd = 1.0E-2), 
               b.hat3=rnorm(1, mean = 0, sd = 1.0E-2),b.hat4=rnorm(1, mean = 0, sd = 1.0E-2),b.hat6=rnorm(1, mean = 0, sd = 1.0E-2),
               taua=rgamma(n = 1, shape = 0.02, rate = 0.02), taub1=rgamma(n = 1, shape = 0.02, rate = 0.02),taub2=rgamma(n = 1, shape = 0.02, rate = 0.02),
               taub3=rgamma(n = 1, shape = 0.02, rate = 0.02),taub4=rgamma(n = 1, shape = 0.02, rate = 0.02))

inits2 <- list(b.hat1=rnorm(1, mean = 0, sd = 1.0E-2), b.hat2=rnorm(1, mean = 0, sd = 1.0E-2), 
               b.hat3=rnorm(1, mean = 0, sd = 1.0E-2),b.hat4=rnorm(1, mean = 0, sd = 1.0E-2),b.hat6=rnorm(1, mean = 0, sd = 1.0E-2),
               taua=rgamma(n = 1, shape = 0.02, rate = 0.02), taub1=rgamma(n = 1, shape = 0.02, rate = 0.02),taub2=rgamma(n = 1, shape = 0.02, rate = 0.02),
               taub3=rgamma(n = 1, shape = 0.02, rate = 0.02),taub4=rgamma(n = 1, shape = 0.02, rate = 0.02))

inits3 <- list(b.hat1=rnorm(1, mean = 0, sd = 1.0E-2), b.hat2=rnorm(1, mean = 0, sd = 1.0E-2), 
               b.hat3=rnorm(1, mean = 0, sd = 1.0E-2),b.hat4=rnorm(1, mean = 0, sd = 1.0E-2),b.hat6=rnorm(1, mean = 0, sd = 1.0E-2),
               taua=rgamma(n = 1, shape = 0.02, rate = 0.02), taub1=rgamma(n = 1, shape = 0.02, rate = 0.02),taub2=rgamma(n = 1, shape = 0.02, rate = 0.02),
               taub3=rgamma(n = 1, shape = 0.02, rate = 0.02),taub4=rgamma(n = 1, shape = 0.02, rate = 0.02))
## Inits
Ngroup<-length(unique(df_coral$stid))
group<-as.factor(df_coral$stid)

Nstid<-3
stid<-rep(seq(1:3),length(unique(df_coral$stid)))

## Data
jd <- list(y=df_coral$cover,Ita.Bleaching=df_coral$ita.bleach.bin,Ita.Nathan.Bleaching=df_coral$ita.nathan.bleach.bin,
           n=nrow(df_coral),ngroup=Ngroup,Group=group,nstid=3,StID=stid)


##Inits

inits1 <- list(b.hat4.0=rnorm(1, mean = 0, sd = 1.0E-2),b.hat5.0=rnorm(1, mean = 0, sd = 1.0E-2),
               taub.0=rgamma(n = 1, shape = 0.02, rate = 0.02))

inits2 <- list(b.hat4.0=rnorm(1, mean = 0, sd = 1.0E-2),b.hat5.0=rnorm(1, mean = 0, sd = 1.0E-2),
               taub.0=rgamma(n = 1, shape = 0.02, rate = 0.02))

inits3 <- list(b.hat4.0=rnorm(1, mean = 0, sd = 1.0E-2),b.hat5.0=rnorm(1, mean = 0, sd = 1.0E-2),
               taub.0=rgamma(n = 1, shape = 0.02, rate = 0.02))

#################### MCMC chains configurations
burnin.values=4000
sample.values= 2000
adapt.values= 1000 
thin.values= 50 

## ! TAKE A LOooooNG TIME TO RUN ! 
mod<-run.jags(model,data= jd, n.chains = 3, method = "parallel",
              inits=list(inits1,inits2,inits3),monitor=c("mu","b4", "b5","b.hat4.0","b.hat5.0","taub.0","taub4","taub5",
                                                         "expected","residual"),burnin = burnin.values, 
              sample = sample.values, adapt = adapt.values,thin= thin.values,
              keep.jags.files = F)

######## Extract model outputs and store them as mcmc.list 

# Model performance
out.perf<- as.mcmc.list(mod, c("expected","residual"))


# Data predictions
out.pred <- as.mcmc.list(mod,var="expected") 

# Parameters estimates
out.beta.hat <-as.mcmc.list(mod,var=c("b.hat1","b.hat2","b.hat3","b.hat4"))
out.beta<- as.mcmc.list(mod, c("b1","b2","b3","b4"))


### Visualisation of MCMC chains convergences

par(mfrow=c(9,3),mar=c(1,1,1,1))
plot(out.beta,density=FALSE,auto.layout = FALSE,ask = F)
dev.off()

par(mfrow=c(9,3),mar=c(1,1,1,1))
plot(out.beta.hat,density=FALSE,auto.layout = FALSE,ask = F)
dev.off()

#gelman.plot(out.beta.hat) --> another MCMC checs=ks

################################################################ Check model performance

Rec.perf<- tidy(out.perf)

df_coral$Expected<-Rec.perf.cumul$Mean[1:nrow(df_coral)]
df_coral$Residuals<-Rec.perf.cumul$Mean[(nrow(df_coral)+1):nrow(Rec.perf)]
df_coral$Nrow<-seq(1:nrow(df_coral_cumul))

# Residuals 
ggplot(df_coral,aes(x=Residuals))+geom_histogram(alpha=0.6)+theme_bw()+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),strip.text = element_text(size=14),
        plot.title = element_text(hjust=0.5,vjust=2,size=15, face="bold"))+
  geom_hline(yintercept = 0,linetype = "longdash")+xlab("Residuals")+ylab("Frequency")

# Predictions vs. residuals

ggplot(df_coral,aes(x=Expected*100,y=Residuals))+geom_point(size=2.8,alpha=0.3)+theme_bw()+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),strip.text = element_text(size=14))+
  geom_hline(yintercept = 0,linetype = "longdash")+xlim(0,60)+xlab("Predicted coral cover (in %)")+ylab("Residuals")

# Observations vs. predictions

ggplot(df_coral,aes(x=cover*100,y=Expected*100))+geom_point(size=2.8,alpha=0.3)+theme_bw()+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),strip.text = element_text(size=14))+
  geom_abline(color="red",linetype="dashed")+xlim(0,60)+ylim(0,60)+xlab("Observed coral cover (in %)")+ylab("Predicted coral cover (in %)")


################################################################ Model goodness-of-fit 

Rec.pred<- tidy(out.pred)

df_coral$Pred<-Rec.pred$Mean
df_coral$Lower<-Rec.pred$X2.5.
df_coral$Upper<-Rec.pred$X97.5.

## Plot coral trajectories with 95% CI from the model and save in your working directory 

df_coral %>% group_by(stid) %>% do({p<-ggplot(data=.) +
  aes(x=year,y=cover*100)  +
  geom_ribbon(aes(ymax=Upper*100, ymin=Lower*100), fill="gray66", alpha=.5)+
  geom_point(size=3.3,color="lightseagreen",alpha=0.6)+
  geom_line(aes(y=Pred*100),size=1)+geom_line(aes(y=cover*100),color="lightseagreen",size=0.4,linetype="3313")+
  theme_bw() +xlab("")+ylab("")+
  ggtitle("Model 1")+ylim(0,50)+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),
        plot.title = element_text(hjust=0.5,vjust=2,size=15, face="bold"))+
  ylim(0,50)+ scale_x_continuous(breaks=c(2012,2013,2014,2015,2016),
                                 labels=c("2012","2013","2014","2015","2016"))
ggsave(p, filename=paste(unique(.$stid),unique(.$reef_name),".pdf"),width=4, height=4, dpi=100,
       path=temp.save)
})

## How many times the model intercerpt the observations? 
df_coral<-df_coral%>%
  mutate(Perf = ifelse(cover >= Lower & cover <=Upper,"In","Out"))

table(df_coral$Perf)/nrow(df_coral)


################################################################ Model parameters 

### Look at the effects of disturbances, at a large spatial scale (Northern Great Barrier Reef)

Rec.beta.hat.cumul<-tidy(out.beta.hat)
param<-c("Ita & Bleaching","Ita & Nathan & Bleaching")

######## Back-transformation from log(odd-ratio) to percentage 

Rec.beta.hat.ind<-Rec.beta.hat.ind%>%mutate(Mean_perc=(exp(Mean)-1)*100)%>%
  mutate(X2.5._perc=(exp(X2.5.)-1)*100)%>%
  mutate(X97.5._perc=(exp(X97.5.)-1)*100)

ggplot(Rec.beta.hat.ind, aes(x=param, y=Mean_perc, fill=param)) + ylab("Effect size (in %)")+ 
  geom_bar(stat="identity",alpha=.6,color="black") +
  geom_errorbar(aes(ymin=X2.5._perc, ymax=X97.5._perc),
                width=.1,position=position_dodge(.9))+theme_bw()+
  theme(axis.text.x = element_text(size=13,angle=30,hjust=1),legend.position="none",
        legend.title = element_text(colour = "black", size = 16, face="bold"), 
        legend.text = element_text(colour = "black", size = 15), 
        panel.grid = element_blank(),plot.title = element_text(hjust=0.5,vjust=2,size=15, face="bold"),
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),axis.title.x=element_text(size=15),
        strip.text=element_text(size=14),strip.background = element_rect(fill="gray98"))+xlab("")+
  geom_hline(yintercept = 0,linetype = "longdash")
