## function to load a model from a filename
# and an annotated file
# do the prediction and export it
library(caret)#,lib.loc='/software/so/el6.3/R-Modules/')
library(randomForest)

predict_me <- function(model_location,file_location,output_name) {
  
  model_to_apply <- readRDS(model_location)
  data_to_predict <-  read.csv(file_location, comment.char = "", check.names = FALSE)
  
  oldnames <- c(names(data_to_predict),'MaxAF','Rank'  ) 
  names(data_to_predict)<-gsub("#", "", names(data_to_predict))
  #set defaults as in the training phase
  data_to_predict$MaxAF<- do.call(pmax,data.frame(data_to_predict$TotalEVSFrequency, data_to_predict$Total1000GenomesFrequency, data_to_predict$ExAC_AF))
  data_to_predict$Cadd2[is.na(data_to_predict$Cadd2) ] <- 0
  data_to_predict$Condel[is.na(data_to_predict$Condel) ] <- 0
  data_to_predict$SegMentDup[is.na(data_to_predict$SegMentDup) ] <- 0
  data_to_predict$PrimatesPhyloP[is.na(data_to_predict$PrimatesPhyloP) ] <- 0
  data_to_predict$PlacentalMammalPhyloP[is.na(data_to_predict$PlacentalMammalPhyloP) ] <- 0
  data_to_predict$PrimatesPhastCons[is.na(data_to_predict$PrimatesPhastCons) ] <- 0
  data_to_predict$PlacentalMammalPhastCons[is.na(data_to_predict$PlacentalMammalPhastCons) ] <- 0
  data_to_predict$Eigen_Phred[is.na(data_to_predict$Eigen_Phred) ] <- 0
  data_to_predict$SIFTScore[is.na(data_to_predict$SIFTScore) ] <- 1
  data_to_predict$MaxAF[is.na(data_to_predict$MaxAF) ] <- 0
  data_to_predict$polyphen2[is.na(data_to_predict$polyphen2) ] <- 0
  data_to_predict$MutAss[is.na(data_to_predict$MutAss) ] <- -5
  data_to_predict$ABB_score[is.na(data_to_predict$ABB_score) ] <- 1
  
  tmp <- predict(model_to_apply,newdata = na.pass(data_to_predict),type='prob')
  data_to_predict$Predicted <- tmp$`1`
  output_data <- data_to_predict[with(data_to_predict, order(-Predicted)), ]
  #save it to CSV
  #output_name <- paste(file_location,'_predicted.csv',sep="")
  names(output_data)<-oldnames

  drops <- c("MaxAF")
  output_data <- output_data[ , !(names(output_data) %in% drops)]


  write.csv(output_data,file=output_name,quote=FALSE,row.names=FALSE)
}
