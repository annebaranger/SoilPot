# Load dataset
library("data.table")
library("tidyr")
library(dplyr)
library(ggplot2)
library(mcp)
library(rstan)
db.clim=fread("output/db_EuForest.csv")

df.traits=read.csv("output/df_trait_filtered.csv") 

pseudoLog10 <- function(x) {asinh(x/2)/log(10)}
#hsm.q01=quantile(db.clim$HSM,0.01,na.rm=TRUE) #because hsm ceilled at -20MPa
db.pres <- db.clim %>%
  rename_with(.cols=everything(),
              tolower) %>% 
  filter(species.binomial%in%c("Pinus sylvestris","Fagus sylvatica","Fraxinus excelsior","Quercus ilex")) %>% 
  #filter(hsm>hsm.q01) %>% 
  filter(!is.na(hsm)) %>% 
  filter(!is.na(fsmwinter)) %>% 
  filter(!is.na(fsmspring)) %>% 
  filter(!is.na(sgdd)) %>% 
  filter(!is.na(wai)) 

hsm.mu=mean(db.pres$hsm,na.rm=TRUE)
hsm.plog.mu=mean(pseudoLog10(db.pres$hsm),na.rm=TRUE)
hsm.sd=sd(db.pres$hsm,na.rm=TRUE)
hsm.plog.sd=sd(pseudoLog10(db.pres$hsm),na.rm=TRUE)
fsm.mu=mean(db.pres$fsmwinter,na.rm=TRUE)
fsm.sd=sd(db.pres$fsmwinter,na.rm=TRUE)

db.pres <- db.pres %>% 
  mutate(hsm=scale(hsm,scale=TRUE,center=TRUE)[, 1],
         hsm.pl10=pseudoLog10(hsm),
         fsmwinter=scale(fsmwinter,scale=TRUE,center=TRUE)[, 1],
         pet=scale(pet,scale=TRUE,center=TRUE)[, 1],
         sgdd=scale(sgdd,scale=TRUE,center=TRUE)[, 1]) 

db.pres.sb <- db.pres %>% 
  # filter(species.binomial%in% c("Fagus sylvatica"))%>% #,"Pinus sylvestris"
  sample_frac(0.01)
# plot explore
# db.pres.sb %>%
  mutate(hsm.ps.sc=pseudoLog10(scale(hsm,scale=TRUE,center=TRUE)),
         hsm.sc.ps=scale(pseudoLog10(hsm),scale=TRUE,center=TRUE)) %>%
  ggplot()+
  geom_point(aes(hsm,hsm.ps.sc),color="green")+
  geom_point(aes(hsm,hsm.sc.ps),color="red")+
#  geom_density(aes(hsm.ps.sc),color="red") + 
  geom_hline(yintercept = -0.07497452,color="green")+
  geom_hline(yintercept = -hsm.plog.mu/hsm.plog.sd,color="red")+
  geom_vline(xintercept=0)

db.pres.sb %>%
  #filter(species.binomial=="Fagus sylvatica") %>% 
  sample_frac(0.1) %>% 
  ggplot(aes(fsmwinter,presence,color=species.binomial))+
  geom_point()+
  geom_vline(xintercept = 0)

# parameters of scaling
hsm.bkpt=-hsm.mu/hsm.sd
fsm.bkpt=-fsm.mu/fsm.sd
hsm.plog.bkpt=pseudoLog10(-hsm.mu/hsm.sd)

# glm fitter
glm.winter=glmmTMB(presence ~  hsm , 
                   data = db.pres.sb,
                   family=binomial(link = "logit")) ## format not working

glm.winter.bis=glm(presence ~  hsm + species.binomial , 
                   data = db.pres.sb,
                   family=binomial(link = "logit"))
# segmented regressions
glm.segmented=segmented(obj=glm.winter.bis,seg.Z=~hsm,npsi=1) 


## models on FSM
# MCP model for one species, one safety margins
library(mcp)
model = list(
  presence ~ fsm,  # changing rate
  ~ 0 
)
prior=list(
  cp_1="dnorm(-0.9199044,1)"
)

data.mcp=data.frame(presence= db.pres.sb$presence,
                    fsm=db.pres.sb$fsmwinter)
fit.mcp.raw= mcp(model, family = bernoulli(),data=data.mcp)
fit.mcp = mcp(model, family = bernoulli(),data=data.mcp,sample = "both", adapt = 10000, prior = prior)
summary(fit.mcp)


#MCP model for several species
model = list(
  presence ~ fsm,  # changing rate
  1+ (1|species.binomial)~ 0 
)
prior=list(
  cp_1="dnorm(-0.9199044,1)"
)

data.mcp=data.frame(presence= db.pres.sb$presence,
                    fsm=db.pres.sb$fsmwinter,
                    species.binomial=as.factor(db.pres.sb$species.binomial))
fit.mcp.raw= mcp(model, family = bernoulli(),data=data.mcp)
fit.mcp = mcp(model, family = bernoulli(),data=data.mcp,sample = "both", adapt = 10000, prior = prior)
summary(fit.mcp)


# stan model, S species, FSM only
data.list <- list(N=dim(db.pres.sb)[1],
                  S=nlevels(factor(db.pres.sb$species.binomial)),
                  presence=db.pres.sb$presence,
                  #hsm=as.numeric(db.pres.sb$hsm),
                  species=as.numeric(factor(db.pres.sb$species.binomial)),
                  fsm=as.numeric(db.pres.sb$fsmwinter),
                  prior_breakpoint=-fsm.mu/fsm.sd,
                  NULL)
fit <- stan(file = "glm_seg_hsm.stan",
            data=data.list,
            iter=1000,
            chains=2)
sum.fit=summary(fit)

##plot results
breakpoint=sum.fit$summary["breakpoint","mean"]
b_fsm=sum.fit$summary["b_fsm","mean"]
plateau=sum.fit$summary["plateau","mean"]
x.fsm=seq(min(data.list$fsm),max(data.list$fsm),len=1000)
p.fsm=inv.logit((x.fsm>breakpoint)*plateau +
  (x.fsm<breakpoint)*(b_fsm * (x.fsm - breakpoint) + plateau))
fit.mean=data.frame(x.fsm,p.fsm)
db.pres.sb %>% 
  ggplot()+
  geom_point(aes(fsmwinter,presence),color="red")+
  geom_line(data=fit.mean,aes(x.fsm,p.fsm))+
  geom_vline(xintercept = data.list$prior_breakpoint)


# stan model, S species, FSM & HSM
data.list <- list(N=dim(db.pres.sb)[1],
                  S=nlevels(factor(db.pres.sb$species.binomial)),
                  presence=db.pres.sb$presence,
                  #hsm=as.numeric(db.pres.sb$hsm),
                  species=as.numeric(factor(db.pres.sb$species.binomial)),
                  fsm=db.pres.sb$fsmwinter,
                  hsm=db.pres.sb$hsm,
                  prior_breakpoint_fsm=-fsm.mu/fsm.sd,
                  prior_breakpoint_hsm=-hsm.mu/hsm.sd,
                  NULL)
fit.2sm <- stan(file = "glm_seg_all.stan",
            data=data.list,
            iter=2000,
            chains=2)
sum.fit=summary(fit.2sm)
save(fit.2sm,file="fit_4sp_all.RData")

##plot results
breakpoint=sum.fit$summary[grepl("breakpoint",rownames(sum.fit$summary)),c("mean","sd")]
b_fsm=sum.fit$summary["b_fsm",c("mean","sd")]
plateau=sum.fit$summary["plateau","mean"]
x.fsm=seq(min(data.list$fsm),max(data.list$fsm),len=1000)
p.fsm=inv.logit((x.fsm>breakpoint)*plateau +
                  (x.fsm<breakpoint)*(b_fsm * (x.fsm - breakpoint) + plateau))
fit.mean=data.frame(x.fsm,p.fsm)
db.pres.sb %>% 
  ggplot()+
  geom_point(aes(fsmwinter,presence),color="red")+
  geom_line(data=fit.mean,aes(x.fsm,p.fsm))+
  geom_vline(xintercept = data.list$prior_breakpoint)



## code GK


x <- seq(-1,1, len = 100)
z <- seq(-1,1, len = 100)
df<- expand.grid(x = x, z = z)
slope1 <- 1

slope2 <- 2

interc <- 0.5
df$y <- (df$x>-0.3)*(interc) + (df$x<=-0.3)*((interc +0.3*slope1)+slope1*df$x) + (df$z>-0.3)*(interc) + (df$z<=-0.3)*((interc+0.3*slope2)+slope2*df$z)

df$y.ilogit <- boot::inv.logit(df$y)

image(x=x,y= z, matrix(df$y, 100, 100))

library(ggplot2)
ggplot(df, aes(x, y, color = z)) + geom_point()
ggplot(df, aes(x, y.ilogit, color = z)) + geom_point()


