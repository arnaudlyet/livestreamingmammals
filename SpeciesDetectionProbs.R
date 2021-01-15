#######################################################################################################################
############  Estimate detection probability per species per Volume - environmental DNA ##################
## Bayesian fremewok
## Arnaud Lyet      
## Started May 20th 2020
## Last update January 15th 2021
#######################################################################################################################
rm(list=ls(all=TRUE))

## packages
library("R2jags")
library("tidyverse")

# model for JAGS
sink("det-prob.txt")
cat("
    
    model{ 
    
    #priors 
    for (j in 1:n.sp){
      p[j] ~ dunif(0,0.2)
    }

    #### model describing detection probability per species per sample
    for(j in 1:n.sp){ 
      for (i in 1:n.samp){ 
        y[i,j] ~ dbern(P[i,j]) # detection event of species i in sample j
        P[i,j] <- 1-((1-p[j])^l[i]) # detection probability p[i] of species i in sample j
      }
    } 
    
    for (j in 1:n.sp){

        P60[j] <- 1-(1-p[j])^60

    }
 
    } #end model description
    
    ",fill = TRUE)
sink()
######################################################


# model for JAGS
sink("det-prob-ct.txt")
cat("
    
    model{ 
    
    #priors 
    for (j in 1:n.sp){
      p[j] ~ dunif(0,0.2)
    }

    #### model describing detection probability per species per sample
    for(j in 1:n.sp){ 
      for (i in 1:n.samp){ 
        y[i,j] ~ dbin(p[j],l[i]) # detection event of species j in sample i
      }
    } 
    
    for (j in 1:n.sp){
 
        P60[j] <- 1-(1-p[j])^60

    }
 
    } #end model description
    
    ",fill = TRUE)
sink()
######################################################



######################################################################
### Data prep for edna and camera traps from edna.full and ct.full
edna.full <- read_csv("edna.sample.csv")
ct.full <- read_csv("ct.sample.csv")

# Create the contingency tables
glimpse(edna.full)

.dna.2018 <- edna.full %>% 
  filter(year==2018) %>%
  mutate(taxa=factor(new_name, levels=list.sp.all)) %>%
  group_by(sample, taxa) %>% 
  summarize(vol=mean(vol), n=1) %>%
  arrange(taxa)

edna.full %>% filter(year==2018) %>% distinct(new_name)

.cont.dna.2018 <- .dna.2018 %>% pivot_wider(names_from = taxa, values_from = n)
.cont.dna.2018[is.na(.cont.dna.2018)] <- 0

n.samp1 <- nrow(.cont.dna.2018)
n.sp1 <- ncol(.cont.dna.2018)-2
y1 <- as.matrix(.cont.dna.2018[3:ncol(.cont.dna.2018)])
l1 <- .cont.dna.2018$vol

.dna.2019 <- edna.full %>% 
  filter(year==2019) %>%
  mutate(taxa=factor(new_name, levels=list.sp.all)) %>%
  group_by(sample, taxa) %>% 
  summarize(vol=mean(vol), n=1) %>%
  arrange(taxa)

.cont.dna.2019 <- .dna.2019 %>% pivot_wider(names_from = taxa, values_from = n)
.cont.dna.2019[is.na(.cont.dna.2019)] <- 0

n.samp2 <- nrow(.cont.dna.2019)
n.sp2 <- ncol(.cont.dna.2019)-2
y2 <- as.matrix(.cont.dna.2019[3:ncol(.cont.dna.2019)])
l2 <- .cont.dna.2019$vol

### camera trap data
glimpse(ct.full)

.ct.2018 <- ct.full %>% 
  filter(year==2018, new_name %in% list.sp.all) %>%
  select(cam, date, effort, new_name) %>% distinct() %>% 
  mutate(taxa=factor(new_name, levels=list.sp.all)) %>%
  group_by(cam, taxa) %>% 
  summarize(effort=mean(effort), n=n()) %>%
  arrange(taxa)

unique(.ct.2018$taxa)

.cont.ct.2018 <- .ct.2018 %>% pivot_wider(names_from = taxa, values_from = n)
.cont.ct.2018[is.na(.cont.ct.2018)] <- 0

n.samp3 <- nrow(.cont.ct.2018)
n.sp3 <- ncol(.cont.ct.2018)-2
y3 <- as.matrix(.cont.ct.2018[3:ncol(.cont.ct.2018)])
l3 <- .cont.ct.2018$effort



.ct.2019 <- ct.full %>% 
  filter(year==2019, new_name %in% list.sp.all) %>%
  select(cam, date, effort, new_name) %>% distinct() %>% 
  mutate(taxa=factor(new_name, levels=list.sp.all)) %>%
  group_by(cam, taxa) %>% 
  summarize(effort=mean(effort), n=n()) %>%
  arrange(taxa)

unique(.ct.2019$taxa)

.cont.ct.2019 <- .ct.2019 %>% pivot_wider(names_from = taxa, values_from = n)
.cont.ct.2019[is.na(.cont.ct.2019)] <- 0

n.samp4 <- nrow(.cont.ct.2019)
n.sp4 <- ncol(.cont.ct.2019)-2
y4 <- as.matrix(.cont.ct.2019[3:ncol(.cont.ct.2019)])
l4 <- .cont.ct.2019$effort



###########################################################################################
#######################    Bayesian analysis for detection    #############################
# MCMC settings
ni <- 11000
nt <- 5
nb <- 1000
nc <- 5

########### edna 2018
n.sp <- n.sp1
n.samp <- n.samp1
y <- y1
l <- l1

data<-list(n.sp=n.sp, n.samp=n.samp, y=y, l=l)
inits<-function(){list(p=runif(n.sp,0,0.01))}

params<-c("P60")



###running the model in rjags
# Call JAGS from R
jagsfit1 <- jagsfit <- jags(data, inits, params, "det-prob.txt", n.chains = nc,
                            n.thin=nt, n.iter=ni, n.burnin=nb, jags.seed = NULL,
                            refresh = (ni-nb)/500, progress.bar = "text")


print(jagsfit,digits=3)


########### VVV edna 2019 VVV ########
n.sp <- n.sp2
n.samp <- n.samp2
y <- y2
l <- l2

data<-list(n.sp=n.sp, n.samp=n.samp, y=y, l=l)
inits<-function(){list(p=runif(n.sp,0,0.01))}

params<-c("P60")



###running the model in rjags
# Call JAGS from R
jagsfit2 <- jagsfit <- jags(data, inits, params, "det-prob.txt", n.chains = nc,
                            n.thin=nt, n.iter=ni, n.burnin=nb, jags.seed = NULL,
                            refresh = (ni-nb)/500, progress.bar = "text")


print(jagsfit,digits=3)


########### Camera trapping 2018
n.sp <- n.sp3
n.samp <- n.samp3
y <- y3
l <- l3

data<-list(n.sp=n.sp, n.samp=n.samp, y=y, l=l)
inits<-function(){list(p=runif(n.sp,0,0.01))}

params<-c("P60")



###running the model in rjags
# Call JAGS from R
jagsfit3 <- jagsfit <- jags(data, inits, params, "det-prob-ct.txt", n.chains = nc,
                            n.thin=nt, n.iter=ni, n.burnin=nb, jags.seed = NULL,
                            refresh = (ni-nb)/500, progress.bar = "text")


print(jagsfit,digits=3)


########### Camera trapping 2019
n.sp <- n.sp4
n.samp <- n.samp4
y <- y4
l <- l4

data<-list(n.sp=n.sp, n.samp=n.samp, y=y, l=l)
inits<-function(){list(p=runif(n.sp,0,0.01))}

params<-c("P60")



###running the model in rjags
# Call JAGS from R
jagsfit4 <- jagsfit <- jags(data, inits, params, "det-prob-ct.txt", n.chains = nc,
                            n.thin=nt, n.iter=ni, n.burnin=nb, jags.seed = NULL,
                            refresh = (ni-nb)/500, progress.bar = "text")


print(jagsfit,digits=3)

# use as.mcmmc to convert rjags object into mcmc.list
#traceplot(jagsfit)



######################################################################
### save results in data frame
res.dna.2018 <- jagsfit1$BUGSoutput$summary
res.dna.2019 <- jagsfit2$BUGSoutput$summary
res.ct.2018 <- jagsfit3$BUGSoutput$summary
res.ct.2019 <- jagsfit4$BUGSoutput$summary

sp.dna.2018<- colnames(.cont.dna.2018)[-c(1:2)]
sp.dna.2019<- colnames(.cont.dna.2019)[-c(1:2)]
sp.ct.2018<- colnames(.cont.ct.2018)[-c(1:2)]
sp.ct.2019<- colnames(.cont.ct.2019)[-c(1:2)]

df.dna.2018 <- data.frame(method="Spatial eDNA",
                           vol=rep("P60", length(sp.dna.2018)), 
                           taxa=sp.dna.2018, res.dna.2018[-nrow(res.dna.2018),])
df.dna.2019 <- data.frame(method="Catchment eDNA",
                          vol=rep("P60",length(sp.dna.2019)), 
                          taxa=sp.dna.2019, res.dna.2019[-nrow(res.dna.2019),])
df.ct.2018 <- data.frame(method="Cam. Trap. 2018",
                          vol=rep("D60", length(sp.ct.2018)), 
                          taxa=sp.ct.2018, res.ct.2018[-nrow(res.ct.2018),])
df.ct.2019 <- data.frame(method="Cam. Trap. 2019",
                          vol=rep("D60", length(sp.ct.2019)), 
                          taxa=sp.ct.2019, res.ct.2019[-nrow(res.ct.2019),])


# get the list of species with the order name

res.tot <- rbind(df.dna.2018,df.dna.2019,df.ct.2018,df.ct.2019) %>%
  select(method, vol, taxa, mean, sd, lwr=X2.5., upr=X97.5.) 

res.tot <- left_join(res.tot,df.sp.gp, by="taxa")


### create a column with detection by which method
.tmp<- res.tot %>% filter(vol=="D60" | vol=="P60") %>% 
  select(method, taxa, mean) %>% pivot_wider(names_from = method, values_from = mean)
.tmp[is.na(.tmp)] <- 0


.tmp <- .tmp %>% mutate(det.dna = `Spatial eDNA`+`Catchment eDNA`, 
         det.ct=`Cam. Trap. 2018` + `Cam. Trap. 2019`,
         det.dna=det.dna/det.dna, det.ct=det.ct/det.ct*3)%>% 
  select(taxa,det.dna,det.ct) %>% replace_na(list(det.dna=0,det.ct=0)) %>%
  mutate(detected=det.dna+det.ct) %>% 
  mutate(.detected=recode(detected, '1'='eDNA', '3'='Camera', '4'='Both'),
         .detected=case_when(taxa=="Canis spp" ~ "Both",
                             taxa=="Weasel spp" ~ "Both",
                             taxa=="Mouse spp" ~ "Both",
                             TRUE ~ .detected
                             )) %>%
  select(taxa, detected=.detected) %>%
  mutate(detected=factor(detected, levels=c("Camera", "Both", "eDNA")))


res.tot <- left_join(res.tot,.tmp, by="taxa")
