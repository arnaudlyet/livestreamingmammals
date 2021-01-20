#########################################################################################################
############  Modelling stream eDNA taxa diversity in relation to    ####################################
############        environmental and sampling covariates            ####################################
## 
## Arnaud Lyet      
## Started October 28th 2020
## Last update December 17th 2020
##
## Code for the GLM shown on Supplementary Figure 3 and Tables 5a and 5b
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
edna.detect <- read_csv(file="edna_detect.csv")


## 4. Preparation of the files for the analysis
# List of species included. Removed Felis sp, Mus musculus and Rangifer tarandus
list.sp.all <- c("Ursus americanus", "Ursus arctos",   
                 "Canis spp", "Canis familiaris", "Canis lupus", "Canis latrans", "Vulpes vulpes",
                 "Lynx canadensis", "Lynx rufus", "Puma concolor", 
                 "Weasel spp", "Mephitis mephitis", "Lontra canadensis", "Gulo gulo", 
                 "Martes americana", "Martes pennanti", "Mustela erminea", "Mustela frenata", 
                 "Mustela nivalis", "Neovison vison",
                 "Alces alces", "Odocoileus hemionus", "Oreamnos americanus", "Ovis canadensis",
                 "Bos taurus", "Equus caballus", 
                 "Lepus americanus", "Ochotona princeps",
                 "Castor canadensis", "Marmota caligata", "Ondatra zibethicus", "Erethizon dorsatum", 
                 "Neotoma cinerea", "Glaucomys sabrinus", "Tamias amoenus", "Tamiasciurus hudsonicus", 
                 "Peromyscus maniculatus", "Myodes gapperi", "Mouse spp", "Arvicolinae spp", 
                 "Sorex cinereus",
                 "Eptesicus fuscus", "Lasionycteris noctivagans", "Myotis spp")


# sum the number of taxa detected per sample
richness.edna <- edna.detect %>%
  filter(taxa %in% list.sp.all) %>%
  group_by(id, rain, aws, vol, time) %>% 
  summarise(n.taxa = sum(det.dna))

df <- richness.edna

# look at the range: 1-26
range(df$n.taxa)

# Definition of the dependent variable. Integer 1 to 26
n.sp.det <- df$n.taxa # nb of taxa detected per sample

# list of covariates
aws <- log(df$aws)
rain <- df$rain
vol <- df$vol
time. <- as.double(df$time/50000)

# parameters to monitor
params <- c("aws", "rain", "vol","time.") 

# creation of the list of models to run
models.list <- mod_lst(params)
models.list <- sapply(models.list, coll_add)

# loop to run the models and store results in model.set
model.set <- list(NULL)
for (t in seq_along(models.list)) {
  .mod <- glm(models.list[t], family=poisson)
  model.set[[t]] <- .mod
}

# Show the model ranked from lower to higher AICc value
head(aictab(model.set, modnames = models.list), 10)

# Beta parameter estimates from the best model
best <- glm(det.dna~aws+rain+vol, family=poisson)
summary(best)


