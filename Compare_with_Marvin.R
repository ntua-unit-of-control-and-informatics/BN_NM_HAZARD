library(bnlearn)
library(yardstick)
library(caret)
library(Rgraphviz)
library(dplyr)

data <- openxlsx::read.xlsx("C:/Users/Vaggelis/Desktop/Vagelis_Bayesian networks for NM Risk Assessment/MarvinsData/Copy of NM_data.xlsx",sheet=1,colNames=TRUE)

#convert all columns to factors
data[] <- lapply(data, as.factor) #lapply returns list so [] keeps the initial data format (data.frame)

#filling in the potential states of variables that were missing from the data set
levels(data$Dissolution)[length(levels(data$Dissolution))+1] <- "25% to 50%"
levels(data$Dissolution)[length(levels(data$Dissolution))+1] <- "50% to 75%"
levels(data$Dissolution)[length(levels(data$Dissolution))+1] <- " >75% "
levels(data$Immunological.effects)[length(levels(data$Immunological.effects))+1] <- "Low"
levels(data$Immunological.effects)[length(levels(data$Immunological.effects))+1] <- "Medium"
levels(data$Immunological.effects)[length(levels(data$Immunological.effects))+1] <- "High"
levels(data$Fibrosis)[length(levels(data$Fibrosis))+1] <- "High"
levels(data$RCNS.effects)[length(levels(data$RCNS.effects))+1] <- "High"
levels(data$Genotoxicity)[length(levels(data$Genotoxicity))+1] <- "High"

str(data)

set.seed(4542352)

N_folds <- 10
N_data <- dim(data)[1]
step <- round(N_data/N_folds)

ug <- bnlearn::empty.graph(names(data))

bnlearn::arcs(ug, check.cycles = TRUE) = matrix(c( "NM.Hazard","Shape", "NM.Hazard","Nanoparticle", "NM.Hazard", "Dissolution",
                                                   
                                                   "NM.Hazard", "Surface.area", "NM.Hazard",
                                                   
                                                   "Surface.charge", "NM.Hazard","Surface.coatings", "NM.Hazard", "Immunological.effects", "NM.Hazard",  
                                                   
                                                   "Surface.reactivity", "NM.Hazard", "Aggregation", "NM.Hazard", 
                                                   
                                                   "Particle.size", "NM.Hazard","Administration.route", "NM.Hazard",
                                                   
                                                   "Study.type", "NM.Hazard","Cytotoxicity", "NM.Hazard",
                                                   
                                                   "Neurological.effects", "NM.Hazard", "Pulmonary.effects", "NM.Hazard",
                                                   
                                                   "Fibrosis", "NM.Hazard", "RCNS.effects", "NM.Hazard", "Genotoxicity", "NM.Hazard",
                                                   
                                                   "Inflammation", "Nanoparticle","Shape", "Nanoparticle", "Dissolution", "Nanoparticle", "Immunological.effects",
                                                   
                                                   "Nanoparticle", "Surface.reactivity", "Nanoparticle", "Neurological.effects",
                                                   
                                                   "Nanoparticle", "Surface.coatings","Nanoparticle", "Surface.charge",
                                                   
                                                   "Nanoparticle", "Administration.route", "Nanoparticle", "Fibrosis",
                                                   
                                                   "Shape","Genotoxicity","Surface.area", "Neurological.effects",
                                                   
                                                   "Surface.coatings", "Surface.area", "Surface.coatings", "Particle.size", 
                                                   
                                                   "Surface.coatings", "Cytotoxicity", "Surface.coatings", "Pulmonary.effects", 
                                                   
                                                   "Surface.coatings", "Aggregation", "Surface.coatings", "Study.type",
                                                   
                                                   "Pulmonary.effects", "Inflammation","Inflammation","RCNS.effects"),
                                                
                                                ncol = 2, byrow = TRUE, dimnames = list(c(), c("from", "to")))



### Structural_EM with Prior Structure and Marvin's data set split

train <- data
test <- openxlsx::read.xlsx("C:/Users/Vaggelis/Desktop/Vagelis_Bayesian networks for NM Risk Assessment/MarvinsData/Copy of NM_data.xlsx",sheet=3,colNames=TRUE)
test[] <- lapply(test, as.factor)

for (j in 1:dim(test)[2]){
    levels(test[,j]) <- levels(train[,j])     
}


  learned <-bnlearn::structural.em(train, maximize = "tabu", maximize.args=list(score= "aic"), 
                                   start=ug, fit = "bayes",   impute = "bayes-lw", impute.args = list(n=2000), return.all = TRUE)
  
  g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(learned$fitted))
  graph::nodeRenderInfo(g) <- list(fontsize=30)
  Rgraphviz::renderGraph(g)
  
  training.fit <- learned$fitted
  
  N_test <- dim(test)[1]
  pred <- rep(0, N_test)
  pred <- as.factor(pred) 
  
  observed.df <- openxlsx::read.xlsx("C:/Users/Vaggelis/Desktop/Vagelis_Bayesian networks for NM Risk Assessment/MarvinsData/Copy of NM_data.xlsx",sheet=4,colNames=TRUE)
  observed <- as.factor(observed.df[,1])
  
  levels(observed) <- c("High",   "Medium", "Low",   "None", "NA")
  levels(pred) <- levels(observed)

for (j in 1:N_test){
 
    pred[j] <- as.character(predict(training.fit, data=test[j,!is.na(test[j,])], node= "NM.Hazard",prob = T, method = "bayes-lw", n=40000))
    
    if(is.na(pred[j])){
      pred[j] <- "NA"
    }
    
}

acc <- sum(pred==observed)/N_test 
mcc <- yardstick::mcc_vec(observed, pred)

df <- data.frame(obs = observed, pred = pred)
confmatList<- df %>%yardstick::conf_mat(obs, pred)

print(acc)

autoplot(confmatList, type = "heatmap")+
scale_fill_gradient(low="#D6EAF8",high = "#2E86C1")+
theme(legend.position = "right")


