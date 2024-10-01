###The script contains the code to train and test the sample-specific XGBoost model, predict combination responses and calculate the t-NSE scores.
###load the necessary packages

pkgs <- c("dplyr","Seurat","HGNChelper", "copykat","readr","ggplot2", "parallel", "HGNChelper", "GSVA", "xgboost", "caret", "ModelMetrics")
lapply(pkgs, library, character.only = T)


###replace the working directory replace '/path/to/working/directory/' with the desired path
path_to_working_directory <- '/path/to/working/directory/'
setwd(dir = path_to_working_directory)

download.file(url = 'https://github.com/yingjchen/RR-AML/archive/refs/heads/main.zip', destfile = 'RR-AML-main.zip')
unzip('RR-AML-main.zip') 
setwd(dir = file.path(path_to_working_directory, 'RR-AML-main'))


##### Step 1: load and process the scRNA-seq data #####
###normalized single cell data with cell types annotated with ScType, and defined malignant and non-malignnat cells
###Take the relapsed AML2 sample as an example 
Expression_data <- readRDS( './exampleData/scAML2R_re.rds' )    ##from zenodo ??
###UMAP showing cell type identification with scType (https://github.com/IanevskiAleksandr/sc-type)
DimPlot(Expression_data, reduction = "umap", label = !0, repel = !0, group.by = 'customclassif') +
  xlab('UMAP1') + ylab('UMAP2') + theme_classic()
ggsave('./Figures/umap_sctype.png',  width = 10, height = 10, dpi = 300)

###check gene symbols
scale_data <- Expression_data@assays[["RNA"]]@scale.data
dim(scale_data)
g_symb <- checkGeneSymbols(rownames(scale_data)) 
scale_data <- scale_data[!is.na(g_symb$Suggested.Symbol), ]; g_symb <- g_symb[!is.na(g_symb$Suggested.Symbol), ]
rownames(scale_data) <- g_symb$Suggested.Symbol; rownames_scale_data <- rownames(scale_data)

##### Optional CopyKAT analysis #####
###identify aneuploid/diploid cells using CopyKAT tool (https://github.com/navinlabcode/copykat)
basecells = c( "Naive CD8+ T cells","Naive CD4+ T cells","Memory CD8+ T cells","Memory CD4+ T cells","Effector CD8+ T cells","Effector CD4+ T cells", "CD4+ NKT-like cells", "CD8+ NKT-like cells")
exp.rawdata <- as.matrix(Expression_data@assays$RNA@counts) # generating this UMI count matrix from 10X output.
normcellnames <- as.character(Cells(Expression_data)[Expression_data@meta.data$customclassif %in% basecells])
copykat.test <- copykat(rawmat=exp.rawdata, id.type ="S", ngene.chr = 5, win.size = 25, KS.cut = 0.1, sam.name = "test", distance = "euclidean", norm.cell.names = normcellnames,output.seg = "FLASE", plot.genes = "TRUE", genome = "hg20",n.cores = 1)

pred.test <- read_delim("./test_copykat_prediction.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
CNA.test <- read_delim("./test_copykat_CNA_results.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnameCNA = colnames(CNA.test)
for(i in 4:ncol(CNA.test)){
  tmp = strsplit(colnames(CNA.test)[i], ".1")
  colnameCNA[i] = paste0(tmp, '-1')
}
length(colnameCNA)
colnames(CNA.test) = colnameCNA

pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]

Expression_data@meta.data$copykatclass = NA
#UMAP showing aneuploid and dipoid cell classification based on CopyKAT
for(j in 1:length(Cells(Expression_data))){
  if (Cells(Expression_data)[j] %in% pred.test$cell.names){
    Expression_data@meta.data$copykatclass[j] = as.character(pred.test$copykat.pred[pred.test$cell.names == Cells(Expression_data)[j]])
  }
}

print('--------ploting--------')

c1_aneuploid = Cells(Expression_data)[Expression_data@meta.data$copykatclass %in% 'aneuploid']
c1_diploid= Cells(Expression_data)[Expression_data@meta.data$copykatclass == 'diploid']


DimPlot(Expression_data, reduction = "umap",label = F, repel = TRUE,  group.by = 'customclassif',
        cells.highlight = list(c1_aneuploid,c1_diploid), pt.size = 1, sizes.highlight = 1,, 
        cols.highlight = c('#08bf2f','#4C78B9'), label.size = 0) +
  guides(color = guide_legend(override.aes = list(size = 6), ncol =1)) +
  scale_color_manual(labels = c('Aneuploid', 'Diploid', NA), values = c('#4C78B9','#08bf2f', '#bcbcbc'), breaks = c('Group_1', 'Group_2', 'Unselected')) + 
  ggtitle( '') +  xlab('UMAP1') + ylab('UMAP2') + labs(colour = "") + theme_classic()
ggsave( './Figures/umap_copykat.png',  width = 10, height = 10, dpi = 300)
##### Optional CopyKAT analysis #####


##### Step 2: process drug-target interactions #####
###load the compound information, including the drug sensitivity scores (DSS) and drug targets 
path_to_DrugInfo <-  './exampleData/exampleData_DrugInfo_AML2R.csv' 
dss_aml1 <- read.csv(path_to_DrugInfo, header = T,sep = ',', check.names = F)

###remove drugs without target information
dss_aml1 <- dss_aml1[dss_aml1[["TargetSet"]]!="",]; dim(dss_aml1)

###check gene symbols
gs_ <- lapply(dss_aml1$Drug, function(d_){ 
  genesInSelectedSets = as.character(dss_aml1[dss_aml1$Drug == d_, "TargetSet"])
  
  # Find gene synonyms
  geneSynonymsMix = unique(na.omit(checkGeneSymbols(unique(unlist(strsplit(genesInSelectedSets,","))))$Suggested.Symbol))
  
  # Genes to keep
  GeneIndToKeep = rownames_scale_data %in% geneSynonymsMix; #keep the genes of scRNA data which can be found in the targeted genes from DSRT data
  ids_ = rownames_scale_data[GeneIndToKeep]; ids_
})


dss_aml1 = dss_aml1[!sapply(gs_, function(i) length(i)==0),]
gs_ = gs_[!sapply(gs_, function(i) length(i)==0)]
names(gs_) = as.character(dss_aml1$Drug)


##### Step 3: prepare the binarized drug/combo-gene interaction matrix; generate expression enrichment matrices for compounds and combinations #####
###single compounds###
gene_union <- unique(unlist(gs_))

###initialize an empty matrix
drug_gene_matrix <- matrix( 0, nrow = length(gs_), ncol = length(gene_union), dimnames = list(names(gs_), gene_union) )

##populate the matrix
for (i in 1:length(gs_) ){
  drug_gene_matrix[i, gs_[[i]]] <- 1
}

##generate expression enrichment matrix
drug_cell_enrichMat <- drug_gene_matrix %*% scale_data[gene_union, ]
###single compounds###


###combinations###
AllCombinations = expand.grid(names(gs_),names(gs_),stringsAsFactors = FALSE)
AllCombinations = AllCombinations[AllCombinations$Var1 != AllCombinations$Var2, ]


AllCombinations$merged = sapply(1:nrow(AllCombinations), function(i) paste0(sort(unlist(AllCombinations[i,])), collapse = ","))  #131406
AllCombinations$Var1 = AllCombinations$Var2 = NULL; AllCombinations_unique = unique(AllCombinations$merged)  #remove the replicated combos 65703

###target sets of combinations
gs_combis = sapply(1:length(AllCombinations_unique), function(i_){
  compounds_s_ = unlist(strsplit(AllCombinations_unique[i_], ","))
  unique(unname(unlist(c(gs_[names(gs_) == compounds_s_[1]], gs_[names(gs_) == compounds_s_[2]]))))
})

names(gs_combis) = AllCombinations_unique

##initialize an empty matrix
combo_gene_matrix <- matrix( 0, nrow = length(gs_combis), ncol = length(gene_union), dimnames = list(names(gs_combis), gene_union) )

##populate the matrix
for (i in 1:length(gs_combis) ){
  combo_gene_matrix[i, gs_combis[[i]]] <- 1
}

combo_cell_enrichMat <- combo_gene_matrix %*% scale_data[gene_union, ]
dim(combo_cell_enrichMat) 
###combinations###


##### Step 4: monotherapy response prediction #####
processed_data = as.data.frame(drug_cell_enrichMat); processed_data$labeloutput = as.numeric(dss_aml1$dss) 

###use grid search and 5-fold cross validation to fine-tune the xgboost model
###parameter set
des <- expand.grid(
  colsample_bytree = seq(0.3, .8, length.out = 5), 
  subsample = seq(0.5, 1.0, length.out = 5), 
  eta = seq(0.01, 0.3, length.out = 5),
  min_child = 1,
  max_depth = 6,
  lambda = 1
  #min_child = seq(3, 4, by = 1),
  #max_depth = seq(3, 4, by = 1),
  #lambda = c(1, 2). #(1,2)
  #nrounds = 2**seq(6, 10, by = 1)
)

CORvalgl <<- list()

###Source code: https://github.com/IanevskiAleksandr/scComb
###generate an objective function and start training
obj.fun = function(x) {
  
  colsample_bytree = x[1]; subsample = x[2]; eta = x[3]; minchild = x[4]; maxdepth = x[5]; lambda = x[6]; MAE_ <- 0; RMSE_ <- 0; 
  COR_ <- 0; COR2_ <- 0; CORval = as.data.frame(cbind(rep(NA, nrow(processed_data)),rep(NA, nrow(processed_data)),rep(NA, nrow(processed_data))))
  
  # repeated CV
  for(repCv in 1:3){
    
    flds1 <- caret::createFolds(1:nrow(processed_data), k = 5, list = T, returnTrain = F);
    
    for(k in 1:length(flds1)){
      testData <- processed_data[flds1[[k]],]; trainData <- processed_data[-flds1[[k]],];
      
      # set N cores = detectCores()-1, feel free to change
      fit = xgboost(data=data.matrix(trainData[, -which(names(trainData) == "labeloutput")]),label = trainData$labeloutput, verbose = F, 
                    nthread = (parallel::detectCores()-1), nrounds = 1024, 
                    params=list(objective = "reg:squarederror", max.depth = maxdepth[[1]], eta=eta[[1]], lambda = lambda[[1]], 
                                min_child_weight = minchild[[1]],
                                colsample_bytree = colsample_bytree[[1]], subsample = subsample[[1]]))
      
      ypred = predict(fit, data.matrix(testData[, -which(names(testData) == "labeloutput")]));
      
      MAE_ <- MAE_ + ModelMetrics::mae(actual = testData$labeloutput, predicted = ypred);
      RMSE_ <- RMSE_ + ModelMetrics::rmse(actual = testData$labeloutput, predicted = ypred);
      COR_ <- COR_ + cor(ypred,testData$labeloutput)
      COR2_ <- COR2_ + cor(ypred,testData$labeloutput, method = "spearman")
      CORval[flds1[[k]],repCv] = ypred
    }
  }
  
  CORvalgl <<- append(CORvalgl,list(CORval))
  paste0(c(MAE_, RMSE_, COR_, COR2_), collapse = ",")
}


###tune the model based on single drugs
start.time <- Sys.time()
gc(T); des$y = apply(des, 1, obj.fun)
message("Finnished Training and fine-tuning in ", round(Sys.time() - start.time, 2), units(Sys.time() - start.time))

###extract all metrics
cols_des = c()
for(i in 1:length(des$y)){ cols_des = rbind(cols_des, as.numeric(as.character(unlist(strsplit(des$y[i],","))))) }
cols_des = as.data.frame(cols_des); colnames(cols_des) = c("MAE","RMSE","PCC","SCC")
des = cbind(des, cols_des)


###pick best 3 models based on pearson correlation coefficient (PCC)
models = des[order(des$PCC),][1:3,]
CORvalgl_top = CORvalgl[order(des$PCC)][1:3]

###get predicted values (inside each cross-validation run, i.e. out-of-fold predictions)
pred_CV = do.call("cbind",lapply(CORvalgl_top, function(i){
  Biobase::rowMedians(data.matrix(i))
}))

# just check how well average out-of-fold predictions correlate with real drug responses
plot(rowMeans(pred_CV), as.numeric(as.character(processed_data$labeloutput)))
ggsave( './Figures/scatter_AML2R_withoutCP.png',  width = 10, height = 10, dpi = 300)
cor(rowMeans(pred_CV), as.numeric(as.character(processed_data$labeloutput)))



##### Step 5: conformal prediction #####
###https://github.com/isidroc/conformal?tab=readme-ov-file
###get absolute errors between average out-of-fold predictions and real response
err_ = abs(Biobase::rowMedians(pred_CV) - as.numeric(as.character(processed_data$labeloutput)))
processed_data_err = processed_data; processed_data_err$labeloutput = err_;

CORvalgl_err <<- list()

###parameter set
des <- expand.grid(
  colsample_bytree = seq(0.3, .8, length.out = 5), 
  subsample = seq(0.5, 1.0, length.out = 5), 
  eta = seq(0.01, 0.3, length.out = 5),
  min_child = 1,
  max_depth = 6,
  lambda = 1
  #min_child = seq(3, 4, by = 1),
  #max_depth = seq(3, 4, by = 1),
  #lambda = c(1, 2). #(1,2)
  #nrounds = 2**seq(6, 10, by = 1)
)

###generate an objective function and start training
obj.fun = function(x) {
  
  colsample_bytree = x[1]; subsample = x[2]; eta = x[3]; minchild = x[4]; maxdepth = x[5]; lambda = x[6]; RMSE_ <- 0; 
  CORval_err = as.data.frame(cbind(rep(NA, nrow(processed_data_err)),rep(NA, nrow(processed_data_err)),rep(NA, nrow(processed_data_err))))
  
  # repeated CV
  for(repCv in 1:3){
    
    flds1 <- caret::createFolds(1:nrow(processed_data_err), k = 5, list = T, returnTrain = F);
    
    for(k in 1:length(flds1)){
      testData <- processed_data_err[flds1[[k]],]; trainData <- processed_data_err[-flds1[[k]],];
      
      # set N cores = detectCores()-1, feel free to change
      fit = xgboost(data=data.matrix(trainData[, -which(names(trainData) == "labeloutput")]),label = trainData$labeloutput, verbose = F, 
                    nthread = (parallel::detectCores()-1), nrounds = 1024, 
                    params=list(objective = "reg:squarederror", max.depth = maxdepth[[1]], eta=eta[[1]], lambda = lambda[[1]], 
                                min_child_weight = minchild[[1]],
                                colsample_bytree = colsample_bytree[[1]], subsample = subsample[[1]]))
      
      ypred = predict(fit, data.matrix(testData[, -which(names(testData) == "labeloutput")]));
      
      RMSE_ <- RMSE_ + ModelMetrics::rmse(actual = testData$labeloutput, predicted = ypred);
      
      CORval_err[flds1[[k]],repCv] = ypred
    }
  }
  
  CORvalgl_err <<- append(CORvalgl_err, list(CORval_err))
  RMSE_ # since we want to focus on large outliers
}



###tune the model based on single drugs
gc(T); des_err$y = apply(des_err, 1, obj.fun)

models_err = des_err[order(des_err$y),][1,]
CORvalgl_top_err = CORvalgl_err[order(des_err$y)][1]
pred_CV_err = do.call("cbind",lapply(CORvalgl_top_err, function(i){
  Biobase::rowMedians(data.matrix(i))
}))

###calculate alpha (conformity scores)
confidence_level = 0.8 # confidence level used for conformity score (0.8 was used as a threshold in the paper), 
alpha <- abs(processed_data$labeloutput - pred_CV[,1]) / pred_CV_err[,1]
alphas <- (sort(alpha))

###predictions to remove (with the largest conformity scores, badly captured by error model)
pred_to_remove <- grep(paste(rownames(processed_data)[alpha > quantile(alpha, confidence_level)],collapse="|"), 
                       rownames(combo_cell_enrichMat), value= !0)


wh = alpha < alphas[length(alphas)*confidence_level]
plot(rowMeans(pred_CV)[wh], as.numeric(as.character(processed_data$labeloutput))[wh])
ggsave( './Figures/scatter_AML2R_withCP.png',  width = 10, height = 10, dpi = 300)
cor(rowMeans(pred_CV)[wh], as.numeric(as.character(processed_data$labeloutput))[wh])
###Conformal end###


##### Step 6: combination response prediction #####
###train three models
models_ <- lapply(1:3, function(i){
  fit_ = xgboost(data=data.matrix(processed_data[, -which(names(processed_data) == "labeloutput")]),label = processed_data$labeloutput, verbose = F, nrounds=1024, 
                 nthread = (parallel::detectCores()-1), params=list(objective = "reg:squarederror", max.depth=models_err[1,"max_depth"], eta=models_err[1,"eta"], 
                                                                    min_child_weight = models_err[1, "min_child"], lambda = models_err[1,"lambda"], subsample = models_err[1,"subsample"],
                                                                    colsample_bytree = models_err[1,"colsample_bytree"]))
  fit_
})
gc(T)

###top features (cells) used by the first model
importance_matrix <- xgb.importance(colnames(processed_data[, -which(names(processed_data) == "labeloutput")]), model = models_[[1]])
xgb.plot.importance(importance_matrix, rel_to_first = TRUE, xlab = "Relative importance")

###UMAP showing the feature importance for  each cell   ---> to be modified 
#Expression_data@meta.data$featureImportance <- c(0)
# Expression_data@meta.data$featureImportance <- NULL  
# Expression_data@meta.data$featureImportance[Cells(Expression_data) %in% importance_matrix$Feature] <- importance_matrix$Importance
# FeaturePlot(Expression_data, features = 'featureImportance', pt.size = 3, order = T) +
#   xlab('UMAP1') + ylab('UMAP2') + theme_classic()
# ggsave(paste0(path0, './Figures/umap_featureImportance.png'),  width = 10, height = 10, dpi = 300)

###predict best combinations
ypredComb = sapply(1:length(models_), function(i){
  predict(models_[[i]], combo_cell_enrichMat)
})  
Combinations_ = data.frame(pred = rowMeans(ypredComb), combis = rownames(combo_cell_enrichMat), stringsAsFactors = !1)
rownames(Combinations_) <- Combinations_$combis

###remove low quality combination predictions (based on conformal predictions)
Combinations_ = Combinations_[!(rownames(Combinations_) %in% pred_to_remove),]
combo_cell_enrichMat = combo_cell_enrichMat[!(rownames(combo_cell_enrichMat) %in% pred_to_remove),]


##### Step 7: combination selection: HSA calculation #####
###Calculate HSA synergy score for combinations
Combinations_$HSA_exp = sapply(strsplit(Combinations_$combis, "\\,"), function(i){
  max(as.numeric(dss_aml1[dss_aml1$drug == i[[1]], "dss"]), as.numeric(dss_aml1[dss_aml1$drug == i[[2]], "dss"]))
})
Combinations_$Synergy = Combinations_$pred - Combinations_$HSA_exp


##### Step 8: combination selection: t-NSE calculation #####
###selective toxicity (estimation), we want to take combinations with the largest t-NSE difference between the malignant cell and non-malignant cells
###based on combo_cell_enrichMat
combo_cell_enrichMat = combo_cell_enrichMat[rownames(Combinations_), ]
gss <- as.data.frame(combo_cell_enrichMat)
###MinMax scaler
###For each sample, perform Minmax scaler on expression enrichment scores of each combo over all the cells 
gss_scale <- as.data.frame(t(apply(gss, 1 ,\(x) (x-min(x))/(max(x)-min(x)))))

###t-NSE score for each cell type
tnse <- data.frame(row.names = rownames(gss))

for (i in unique(Expression_data$customclassif)){
  tnse[, i] = rowSums(as.matrix(gss_scale[, Expression_data$customclassif == i]))/sqrt(sum(Expression_data$customclassif == i))
}


###t-NSE score for cancer and normal cells
cancer_celltype <- unique(Expression_data@meta.data[["customclassif"]][Expression_data@meta.data[["malignant"]] ] )
normal_celltype <- unique(Expression_data@meta.data[["customclassif"]][Expression_data@meta.data[["nonmalignant"]] ] )
tnse[, 'CancerClusters'] = rowMeans(as.matrix(tnse[, cancer_celltype ]))
tnse[, 'NormalClusters'] = rowMeans(as.matrix(tnse[, normal_celltype ]))
Combinations_[,"selective_toxicity"] <- tnse$CancerClusters - tnse$NormalClusters


##### Step 9: final selection of combinations for validation #####
###select combinations with expected synergy HSA > 5
Combinations_ = Combinations_[Combinations_$Synergy > 5]
###select combinations with highest expected efficacy
Combinations_ = Combinations_[Combinations_$pred > quantile(Combinations_$pred, .9),]
###select combinations with higher t-NSE score differences between the cancer and normal cells 
Combinations_ = Combinations_[Combinations_$selective_toxicity > quantile(Combinations_$selective_toxicity, .5),]