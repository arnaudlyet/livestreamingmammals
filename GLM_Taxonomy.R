#######################################################################################################################
############  Look at the spatial relation between camera trap and eDNA ##################
## 
## Arnaud Lyet      
## Started October 28th 2020
## Last update December 17th 2020
##
## Code for the GLM shown on Figure 4 and Supplementary Tables 2a and 2b
## 1. Packages
## 2. Function to create list of all possible combinations from a list of parameters
## 3. Loading the row edna and camera trapping tables
## 4. Preparation of the files for the analysis
## 5. Modeling and model selection
#######################################################################################################################
rm(list=ls(all=TRUE))


## 1. packages
library("tidyverse")
library("modelr")
library("AICcmodavg")


## 2. function to create list of all possible combinations from a list of parameters
mod_lst <-  function(params) {
  models <- (NULL)
  for(l in seq_along(params)) {
    .models <- combn(params, m=l, simplify=FALSE)
    models <- c(models, .models)
  }
  c(1,models)
} 

# functions to create the models in additive or interactive manner
coll_add <- function(x) {paste("det.dna~",paste(x, collapse="+"), sep="")}
coll_int <- function(x) {paste("det.dna~",paste(x, collapse="*"), sep="")}

### inverse logit
invlogit <- function (x) {1 / (1 + exp(-x))}

### logit
logit <- function (x) {log(x/(1-x))}

#### END functions


## 3. Loading the files containing the eDNA and camera trappind data
edna.ct.match <- read_csv(file="edna_ct_match.csv")
edna.detect <- read_csv(file="edna_detect.csv")
ct.2018 <- read_csv(file="ct_2018.csv")
ct.2019 <- read_csv(file="ct_2019.csv")


## 4. Preparation of the files for the analysis
# Matching eDNA sites/samples with camera traps species detections within catchment
# Spatial design vs CT 2018
edna.ct.match.2018 <- edna.ct.match %>% filter(year==2018)
output.2018 <- NULL

for (i in 1:nrow(edna.ct.match.2018)) {
  
  loc <- as.character(edna.ct.match.2018[i,"id"]) # pick code sample eDNA evaluated
  ct.ws <- ct.names[!is.na(as.vector(edna.ct.match.2018[i,5:66]))]
  
  tmp.eff <- ct.2018 %>% 
    filter(cam %in% ct.ws) %>% 
    group_by(cam, effort) %>%  
    summarize()
  
  eff.tot <- sum(tmp.eff$effort)
  
  # temporary tab with species, avg.freq, sample id
  tmp.ws <- ct.2018 %>% 
    filter(cam %in% ct.ws) %>% 
    group_by(new_name2) %>% ### new_name2 instead of species because of Canis 
    summarize(det.ws=n()/eff.tot) %>% # 
    mutate(id=loc, year=2018)
  
  # assemble the final table
  output.2018 <- rbind(output.2018, tmp.ws)
  
}

# Catchment design vs CT 2019
edna.ct.match.2019 <- edna.ct.match %>% filter(year==2019)
output.2019 <- NULL

for (i in 1:nrow(edna.ct.match.2019)) {
  
  loc <- as.character(edna.ct.match.2019[i,"id"]) # pick code sample eDNA evaluated
  ct.ws <- ct.names[!is.na(as.vector(edna.ct.match.2019[i,5:66]))]
  
  tmp.eff <- ct.2019 %>% 
    filter(cam %in% ct.ws) %>% 
    group_by(cam, effort) %>%  
    summarize()
  
  eff.tot <- sum(tmp.eff$effort)
  
  # temporary tab with species, avg.freq, sample id
  tmp.ws <- ct.2019 %>% 
    filter(cam %in% ct.ws) %>% 
    group_by(new_name2) %>% ### new_name2 instead of species because of Canis 
    summarize(det.ws=n()/eff.tot) %>% # 
    mutate(id=loc, year=2019)
  
  # assemble the final table
  output.2019 <- rbind(output.2019, tmp.ws)
  
}

# merging tables 2018 and 2019
output.ct.det <- rbind(output.2018, output.2019) %>%
  select(id, taxa=new_name2, det.ws)



# Selection of the taxa to include in the analysis and taxonomic group assignment

# List of species included: removed domestic species, grouped canis, grouped small ground rodents
list.sp.ct <- c("Ursus americanus", "Ursus arctos",   
                "Canis spp", "Vulpes vulpes", # all species in the genus Canis under Canis spp 
                "Lynx canadensis", "Lynx rufus", "Puma concolor", 
                "Mephitis mephitis", "Gulo gulo", "Martes americana",  
                "Martes pennanti", "Mustela erminea", "Mustela frenata" ,"Mustela nivalis",
                "Alces alces", "Odocoileus hemionus", # domestic species removed c("Bos taurus", "Equus caballus") 
                "Lepus americanus", "Marmota caligata", "Erethizon dorsatum", "Neotoma cinerea", 
                "Tamias amoenus", "Mouse spp", # c("Peromyscus maniculatus","Myodes gapperi","Arvicolinae spp") grouped under Mouse spp
                "Glaucomys sabrinus", "Tamiasciurus hudsonicus") 

df.dna <- edna.detect %>%  
  mutate(taxa = case_when(
    taxa %in% c("Peromyscus maniculatus", 
                "Myodes gapperi", 
                "Arvicolinae spp") ~ "Mouse spp",
    TRUE ~ taxa)) %>% 
  filter(taxa %in% list.sp.ct) %>%
  mutate(group = case_when(
    taxa %in% list.sp.ct[1:2] ~ "Bears",
    taxa %in% list.sp.ct[3:4] ~ "Canids",
    taxa %in% list.sp.ct[5:7] ~ "Felids",
    taxa %in% list.sp.ct[8:14] ~ "Mustelids",
    taxa %in% list.sp.ct[15:16] ~ "Ungulates",
    taxa %in% list.sp.ct[17] ~ "Lagomorphs",
    taxa %in% list.sp.ct[18:22] ~ "Ground rodents",
    taxa %in% list.sp.ct[23:24] ~ "Arboreal rodents")
  )

df.ct <- output.ct.det %>% 
  filter(taxa %in% list.sp.ct)

df <- left_join(df.dna,df.ct, by = c("id","taxa"))
df$det.ws[is.na(df$det.ws)] <- 0

df <- df %>% mutate(group = factor(group, levels=c("Bears","Canids","Felids","Mustelids",
                                                   "Ungulates","Lagomorphs","Ground rodents",
                                                   "Arboreal rodents")))


## 5. Modeling and model selection

# Dependent variable 1/0 format
det.dna <- df$det.dna

# list of covariates
det.ws <- df$det.ws # average detection of species at the catchment
rain <- df$rain # millimetre of rain the day before sampling
aws <- log(df$aws) # log transformed
vol <- df$vol # litre
year <- as.factor(df$year) #factor
gp <- df$group # factor


# parameters to monitor
params <- c("det.ws","rain","aws","vol","gp","year") 

# creation of the list of models to run
models.list <- mod_lst(params)
models.list <- sapply(models.list, coll_add)

# loop to run the models and store results in model.set
model.set <- list(NULL)
for (t in seq_along(models.list)) {
  .mod <- glm(models.list[t], family=binomial)
  model.set[[t]] <- .mod
}

# Show the model ranked from lower to higher AICc value
head(aictab(model.set, modnames = models.list), 10)

# Beta parameter estimates from the best model
best <- glm(det.dna~det.ws+rain+aws+vol+gp, family=binomial)
summary(best)

