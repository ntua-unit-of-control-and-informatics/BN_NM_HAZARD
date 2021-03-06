library(bnlearn)
library(yardstick)
library(caret)
library(Rgraphviz)
library(dplyr)

data<-openxlsx::read.xlsx("C:/Users/Vaggelis/Desktop/Vagelis_Bayesian networks for NM Risk Assessment/MarvinsData/Copy of NM_data.xlsx",sheet=6,colNames=TRUE)

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

set.seed(1298)

N_folds <- 10
N_data <- dim(data)[1]
step <- round(N_data/N_folds)

### Bayesian Network with Constrained based algorithms

acc_best <- rep(0, N_folds)
mcc_best <- rep(0, N_folds)
confmatList <- list()
start <- 1

#10_fold cross validation and prediction for the NM.hazard variable
for(i in 1: N_folds){
  
  select <- start:(start+step-1)
  
  test <- data[select,]
  observed <- test$NM.Hazard
  test[,12:18] <- NA 
  train <- data[-select,]
  
  learned <- bnlearn::gs(train)
  
  learned1 <- cextend(learned, strict = TRUE) 

  g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(learned1))
  graph::nodeRenderInfo(g) <- list(fontsize=30)
  Rgraphviz::renderGraph(g)
  
  training.fit <- bnlearn::bn.fit(learned1, train, method="bayes")
 
  N_test <- dim(test)[1]
  pred <- rep(0, N_test)
  
  pred <- as.factor(pred) 
  observed <- as.factor(observed)
  levels(observed) <- c("High",   "Low",    "Medium", "None", "NA")  
  levels(pred) <- levels(observed)
  
  for (j in 1:N_test){
    if(i == 10 && j >=56){
      break;
    }else{
      pred[j] <- as.character(predict(training.fit, data=test[j,!is.na(test[j,])], 
                                      
                                      node= "NM.Hazard",prob = T, method = "bayes-lw", n=40000))
      if(is.na(pred[j])){
        pred[j] <- "NA"
      }
    }
  }
  if (i != 10){
    acc_best[i] <- sum(pred==observed)/N_test 
    mcc_best[i] <- yardstick::mcc_vec(observed, pred)
    df <- data.frame(obs = observed, pred = pred)
    confmatList[[i]] <- df %>%yardstick::conf_mat(obs, pred)
  }else{
    pred <- pred[1:length(pred)-1]
    observed <- observed[1:length(observed)-1]
    acc_best[i] <- sum(pred==observed)/N_test 
    mcc_best[i] <- yardstick::mcc_vec(observed, pred)
    
    df <- data.frame(obs = observed, pred = pred)
    confmatList[[i]] <- df %>%yardstick::conf_mat(obs, pred)
    
  }
  
  start <- start+step
  
}

#finding the average of Accuracy and MCC
mean_acc_best <- mean(acc_best)
mean_mcc_best <- mean(mcc_best)

#generating the average confusion matrix 
average.confmat<-  confmatList[[1]]
for (i in 1:4){
  for (j in 1:4){
    average.confmat$table[i,j] <- round(mean(c(confmatList[[1]]$table[i,j], confmatList[[2]]$table[i,j], confmatList[[3]]$table[i,j], 
                                               confmatList[[4]]$table[i,j], confmatList[[5]]$table[i,j], confmatList[[6]]$table[i,j],
                                               confmatList[[7]]$table[i,j], confmatList[[8]]$table[i,j], confmatList[[9]]$table[i,j],
                                               confmatList[[10]]$table[i,j] )))
  }
}


autoplot(average.confmat, type = "heatmap")+
  scale_fill_gradient(low="#D6EAF8",high = "#2E86C1")+
  theme(legend.position = "right")
