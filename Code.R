library(caret)
library(infotheo)
library(FSelector)
library(stringr)
library(doParallel)
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
set.seed(1)

filter_features_sets <- function(data, is_cumul, window_size, slice) {
  aoinames = c("Relevant.bars", "Non.relevant.bars", "Text", "Refs", "labels", "Viz", "legend")
  toremove = c("Sc_id", "blinknum", "blinkdurationmax", "blinktimedistancestd", "blinktimedistancemin", 
               "blinktimedistancemean","blinkdurationmean", "blinktimedistancemax", "blinkrate","blinkdurationmin", 
               "blinkdurationtotal", "blinkdurationstd", "length", "numfixations", "numsegments", "doubleclicrate", 
               "sumabspathangles", "sumfixationduration" ,"sumpathdistance" ,"sumrelpathangles", "sumsaccadedistance", 
               "sumsaccadeduration","numevents", "numleftclic", "leftclicrate", "numrightclic", "numdoubleclic", "numsaccades", 
               "numsamples", "numkeypressed", "rightclicrate", "keypressedrate","timetofirstdoubleclic", "timetofirstkeypressed", "timetofirstleftclic",
               "timetofirstrightclic", "timetolastfixation", "Uid")
  
  for(i in 1:length(aoinames)){
    toremove = c(toremove, paste(aoinames, "blinknum", sep="_"))
    toremove = c(toremove, paste(aoinames, "doubleclicrate", sep="_"))
    toremove = c(toremove, paste(aoinames, "numevents", sep="_"))
    toremove = c(toremove, paste(aoinames, "numfixations", sep="_"))
    toremove = c(toremove, paste(aoinames, "numleftclic", sep="_"))
    toremove = c(toremove, paste(aoinames, "numrightclic", sep="_"))
    toremove = c(toremove, paste(aoinames, "numdoubleclic", sep="_"))
    toremove = c(toremove, paste(aoinames, "rightclicrate", sep="_"))
    toremove = c(toremove, paste(aoinames, "timetofirstdoubleclic", sep="_"))
    toremove = c(toremove, paste(aoinames, "timetofirstrightclic", sep="_"))
    toremove = c(toremove, paste(aoinames, "timetofirstleftclic", sep="_"))
    toremove = c(toremove, paste(aoinames, "timetolastdoubleclic", sep="_"))
    toremove = c(toremove, paste(aoinames, "timetolastfixation", sep="_"))
    toremove = c(toremove, paste(aoinames, "timetolastleftclic", sep="_"))
    toremove = c(toremove, paste(aoinames, "timetolastrightclic", sep="_"))
    toremove = c(toremove, paste(aoinames, "totaltimespent", sep="_"))
    toremove = c(toremove, paste(aoinames, "enddistance", sep="_"))
    toremove = c(toremove, paste(aoinames, "endpupilsize", sep="_"))
    toremove = c(toremove, paste(aoinames, "startdistance", sep="_"))
    toremove = c(toremove, paste(aoinames, "startpupilsize", sep="_"))
  }
  data = data[, !(colnames(data) %in% toremove)]
  distance_features = c("enddistance", "maxdistance", "meandistance", "mindistance", "startdistance")
  for (feature in distance_features) {
    data[data[feature] == -1, feature] = 600
    for (aoi in aoinames) {
      feat_name = paste(aoi, feature, sep="_")
      if (feat_name %in% colnames(data)) {
        data[data[feat_name] == -1, feat_name] = 600
      }
    }
  }
  stddev_other_mean_features = c("stddevabspathangles", "stddevdistance", "stddevfixationduration", "stddevpathdistance", "stddevpupilsize", "stddevrelpathangles", 
                                 "stddevsaccadedistance", "stddevsaccadeduration", "stddevsaccadespeed", "meanabspathangles", "meanfixationduration", "meanpathdistance", 
                                 "meanpupilsize", "meanpupilvelocity", "meanrelpathangles", "meansaccadedistance", "meansaccadeduration", "meansaccadespeed")
  for (feature in stddev_other_mean_features) {
    data[data[feature] == -1, feature] = 0
    for (aoi in aoinames) {
      feat_name = paste(aoi, feature, sep="_")
      if (feat_name %in% colnames(data)) {
        data[data[feat_name] == -1, feat_name] = 0
      }
    }
  }  
  if (is_cumul) {
    for (feature in c("timetofirstfixation")) {
      for (aoi in aoinames) {
        feat_name = paste(aoi, feature, sep="_")
        if (feat_name %in% colnames(data)) {
          data[data[feat_name] == -1, feat_name] = window_size * slice
        }
      }
    }
  }
  return(data)
}

fill_results <- function(model_stats, model, results, iter_idx) {
  model_stats[iter_idx, paste(model, "_Accuracy", sep="")] =  results$overall[1]
  model_stats[iter_idx, paste(model, "_Kappa", sep="")] = results$overall[2]
  model_stats[iter_idx, paste(model, "_AccuracyHigh", sep="")] = results$byClass[[1]]
  model_stats[iter_idx, paste(model, "_AccuracyLow", sep="")] = results$byClass[[2]]
  model_stats[iter_idx, paste(model, "_NumLow", sep="")] = sum(testing_data[, dependent[dep]] == 0)
  model_stats[iter_idx, paste(model, "_NumHigh", sep="")] = sum(testing_data[, dependent[dep]] == 1)
  
  return(model_stats)
}

train_test_model <- function(model_name, tune_params, models_stats, feature_imp_fname, fold) {
  fit <- train(fmla,
               data = training_data,
               method = model_name,
               tuneGrid = tune_params,
               #tuneLength = tune_params,
               trControl = ctrl,
               na.action=na.pass,
               metric="Accuracy")
  
  sink(feature_imp_fname, append = TRUE)
  print(paste("Slice ",slice, ", most important features for feature ", dependent[dep], ", run ", numRepeat, ", model ", model_name, sep=""))
  print(varImp(fit))
  sink() 
  res = predict(fit, newdata=testing_data)
  res = confusionMatrix(res, testing_data[,dependent[dep]])
  models_stats = fill_results(models_stats, model_name, res,  fold)
  return (models_stats)
}
weighted_mean_sd <- function(model, pref, models_stats) {
  sumn = sum(models_stats[, paste(model, "_", "Num", pref, sep="")])
  weights = models_stats[, paste(model, "_", "Num", pref, sep="")] / sumn
  wmean = models_stats[, paste(model, "_", "Accuracy", pref, sep="")] %*% weights
  M = sum(weights != 0)
  sdw = sqrt(weights %*% ((models_stats[, paste(model, "_", "Accuracy", pref, sep="")] - wmean)^2) / (sum(weights) * (M-1)/M)) 
  result = paste(result, wmean, sep=",")
  result = paste(result, sdw, sep=",")
  return(result)
}
add_labels_to_features <- function(features_data, avail_users, label_data, labels_msnv_data, user_dependent, msnv_dependent, keep_hard_msnvs) {
  features_data["Uid"] = -1
  curr_usr = -1
  if (cumul_across) {
    for(row in 1:nrow(features_data)) {
      sc_id = toString(features_data[row, "Sc_id"])
      if (endsWith(sc_id, "allsc")) {
        undersc_pos = regexpr("_", sc_id)[1] - 1
        curr_usr = strtoi(substring(sc_id, 1, undersc_pos))
        if (curr_usr %in% avail_users) {
          features_data[row, user_dependent] = label_data[label_data[,"user_id"] == curr_usr, user_dependent]
          features_data[row, "Uid"] = curr_usr
        }
      }
    }
  } else{
    for(row in 1:nrow(features_data)) {
      sc_id = toString(features_data[row, "Sc_id"])
      if (endsWith(sc_id, "allsc")) {
        undersc_pos = regexpr("_", sc_id)[1] - 1
        curr_usr = strtoi(substring(sc_id, 1, undersc_pos))
      } else {
        if (curr_usr %in% avail_users) {
          msnv_id = strtoi(features_data[row, "Sc_id"])
          if (msnv_id == 18 || is.na(msnv_id) || (keep_hard_msnvs && !(msnv_id %in% c(3, 5)))) {
            
          } else {
            features_data[row, user_dependent] = label_data[label_data[,"user_id"] == curr_usr, user_dependent]
            features_data[row, msnv_dependent] = labels_msnv_data[labels_msnv_data[,"user_id"] == curr_usr & labels_msnv_data[,"mmd_id"] == msnv_id, msnv_dependent]
            features_data[row, "Uid"] = curr_usr
          }
        }
      }
    }
  }
  # Clean up the data
  features_data = features_data[features_data[dependent[1]] != -1, ]
  features_data[is.na(features_data)] = -1
  return(features_data)
}

############# SET THESE BEFORE RUNNING #####################
# Folder to save prediction outputs/feature importances
dirout = "C:\\Users\\Alire\\Desktop\\rcode\\results\\" 
# Path to folder with labels
label_dir = "C:\\Users\\Alire\\Desktop\\rcode\\raw_experiment\\Tobii Export\\"
# Path to .txt files containing BarChartLit test results
users_files_path <- "C:\\Users\\Alire\\Desktop\\rcode\\user_data\\first_half\\"
# Features that are independent of specific
user_dependent = c("BarChartLit") # c("BarChartLit", "VerbalWM_longest", "Meara")
# Features dependent on specific MSNVs
msnv_dependent = c()#"mmd_accuracy", "mmd_task_time", "mmd_interest")
# Output of which models to save
models = c("MajorityClass","BstLm")#"rf", "regLogistic", "xgbTree","svmLinear") #LogitBoost,
# Metrics to store
metrics = c("Accuracy","Kappa", "AccuracyHigh", "AccuracyLow")
#############################################################
#==========================================
feat_set = "all"
feat_selection = TRUE
K = 10


# NOTE: if all false, will run cumulative within window
disjoint_window = F # always false
cumul_across = T # T: across task 1 to 15 tasks | F: 29 sec                                         ####### Across/within
full_window = F # always false
keep_hard_msnvs = T # always true

if (disjoint_window) {
  time_windows = c(2000, 5000, 10000)
  window_id = 1
  out = paste(dirout, "predict-taskchar_result_disjoint_window_", time_windows[window_id],".csv", sep="")
  feature_imp_fname = paste(dirout, "disjoint_window.txt",sep="")
} else if (cumul_across) {
  time_windows = c(1000)
  out = paste(dirout, "predict-taskchar_result_accross.csv", sep="")
  feature_imp_fname = paste(dirout, "featimp_cumul_accross.txt",sep="")
} else if (full_window) {
  time_windows = c(1000)
  out = paste(dirout, "predict-taskchar_result_fullwindow.csv", sep="")
  feature_imp_fname = paste(dirout, "featimp_fullwindow.txt",sep="")
} else {
  time_windows = c(1000)
  seqs = c(1, 9, 10, 19, 20, 29, 30, 39, 40, 49, 50, 58)                                          ####### seconds
  machine_id = 1
  left_idx = seqs[2*(machine_id-1)+1]
  right_idx = seqs[2*(machine_id-1)+2]
  feature_imp_fname = paste(dirout, "featimp_cumul_within_",machine_id,".txt",sep="")
  out = paste(dirout, "predict-taskchar_result_cumulative_",left_idx,"_",right_idx,".csv",sep="")
}

sink(out)
header = "Userchar,WindowSize,Slice,Model,Repeat"
for (metric in metrics) {
  header = paste(header, ",", metric, "_", "mean", sep="")
  header = paste(header, ",", metric, "_", "sd", sep="")
}
print(header)
sink()

#==========================================
# Read labels data
labels_filename = paste(label_dir, "user_chars_MMD.csv", sep="")

label_data = read.table(labels_filename, header=T, sep=",")
label_data = na.omit(label_data)


users_files = list.files(users_files_path)
user_barchart_scores =  rep(0, 100)

for (file_name in users_files)  {
  linesplit = strsplit(file_name, "_")
  user_id = strtoi(linesplit[[1]][1])
  conn <- file(paste(users_files_path, file_name, sep=""),open="r")
  lines <-readLines(conn)
  total_score = 0
  num_attempts = 0
  for (i in 2:length(lines)){
    line = lines[i]
    if (str_detect(line, "Score")) {
      num_attempts = num_attempts + 1
      total_score = total_score + strtoi(strsplit(line, ":")[[1]][2])
    }
  }
  user_barchart_scores[user_id] = total_score
  close(conn)
}


labels_msnv_filename = paste(label_dir, "performance_data_msnv_all.csv", sep="")
labels_msnv_data = read.table(labels_msnv_filename, header=T, sep=",")
labels_msnv_data = na.omit(labels_msnv_data)

#Median split stuff
verbal_VM_words_median = median(label_data$VerbalWM_word)
barchart_median = median(user_barchart_scores[user_barchart_scores != 0])

dependent = c(user_dependent, msnv_dependent)

users_in_cumul = c(95, 67, 66, 65, 64, 63, 62, 92, 50, 46, 45, 42, 40, 38, 93, 61, 60, 59, 58, 55, 52, 97, 73, 72, 71, 70, 69, 68, 90, 19, 18, 16, 12,  9,
                   1, 91, 36, 31, 30, 26, 25, 21, 89, 88, 85, 84, 81, 80, 79, 78, 77, 76, 75, 74)
labels_msnv_data = labels_msnv_data[labels_msnv_data$user_id %in% users_in_cumul, ]
for(i in 1:length(msnv_dependent)) {                                                        #Over performance measures
  if(is.numeric(labels_msnv_data[, msnv_dependent[i]])) { #DISCRETIZE LABELS
    dep_median = median(labels_msnv_data[, msnv_dependent[i]])
    for (row in 1:nrow(labels_msnv_data)) {
      numeric_val = labels_msnv_data[row, msnv_dependent[i]]
      if (numeric_val > dep_median) {
        labels_msnv_data[row, msnv_dependent[i]] = 1
      } else {
        labels_msnv_data[row, msnv_dependent[i]] = 0
      }
    }
  }
}

label_data = label_data[label_data$user_id %in% users_in_cumul, ]
for(i in 1:length(user_dependent)) {                                                        #Over cognitive measures
  if(is.numeric(label_data[, user_dependent[i]])) { #DISCRETIZE LABELS
    dep_median = median(label_data[, user_dependent[i]])
    for (row in 1:nrow(label_data)) {
      numeric_val = label_data[row, user_dependent[i]]
      if (numeric_val > dep_median) {
        label_data[row, user_dependent[i]] = 1
      } else {
        label_data[row, user_dependent[i]] = 0
      }
      if (user_dependent[i] == "VerbalWM_longest" && numeric_val == dep_median) {
        if (label_data[row, "VerbalWM_word"] > verbal_VM_words_median) {
          label_data[row, user_dependent[i]] = 1
        }
      } else if (user_dependent[i] == "BarChartLit") {
        uid =  label_data[row, "user_id"]
        numeric_val =  user_barchart_scores[uid]
        if (numeric_val > barchart_median) {
          label_data[row, user_dependent[i]] = 1
        } else {
          label_data[row, user_dependent[i]] = 0
        }
      }
    }
  }
}

# Set of users to keep in dataframe
avail_users = unique(label_data[, "user_id"])

for (window in time_windows) {
  if (cumul_across) { # if across task task_1-task_15
    features_dir = "C:\\Users\\Alire\\Desktop\\rcode\\train_data\\across_tasks_new_refined\\"
    iter_arr = c(3)         # seq(1, 15, 1)                                                      # Across task: # task
  } else if (full_window) {
    features_dir = "C:\\Users\\Alire\\Desktop\\rcode\\train_data\\across_tasks_new_refined\\"
    iter_arr = seq(15, 15, 1)
    features_dir = "C:\\Users\\Alire\\Desktop\\rcode\\train_data\\"
  } else {
    # Cumulative within
    features_dir = "C:\\Users\\Alire\\Desktop\\rcode\\train_data\\cumulative_new_withintask_refined\\"
    iter_arr = seq(left_idx, right_idx, 1)                                                        # within task 1s-29s
  }
  for(slice in iter_arr) {                                                                        # Over each slice
    if (cumul_across || full_window) {
      features_filename = paste(features_dir, "tasks_included_",slice,".tsv", sep="")
    } else {
      features_filename = paste(features_dir, "pruning_", slice*1000,".tsv", sep="")
    }
    
    features_data = read.table(features_filename, header=T, sep="\t")
    # Add label columns
    features_data[dependent] = -1
    features_data = add_labels_to_features(features_data, avail_users, label_data, labels_msnv_data, user_dependent, msnv_dependent, keep_hard_msnvs)
    
    user_ids = unique(features_data$Uid)
    # For all user characteristics
    for(dep in 1:length(dependent))
    {
      # Keep only 1 dependent
      dependent2 = dependent[dependent != dependent[dep]]
      features_data2 = features_data[, !(colnames(features_data) %in% dependent2)]
      models_stats = data.frame()
      # Statistics
      table_dependent = table(features_data[dependent[dep]])
      maj_value = names(table_dependent[table_dependent == max(table_dependent)])
      maj_acc = nrow(features_data[features_data[,dependent[dep]] == maj_value,]) / nrow(features_data)
      
      #Repeated CV--------------------------------------------------------------------------------------------
      # Repeat fold splits multiple times to avoid lucky splits
      num_repeats = 1
      K = 53
      outerloop_train_sets = list()
      outerloop_test_sets = list()
      Sum = 0
      
      # print(slice)
      
      
      for (numRepeat in 1:num_repeats) {
        user_ids = user_ids[sample(length(user_ids))]
        folds = cut(seq(1, length(user_ids)), breaks=K, labels=F)
        
        # Go over each fold
        for(fold in 1:K) {
          
        print(fold)

          models_stats[fold, "MajorityClass_Accuracy"] = maj_acc
          
          # Extract train and validation data
          testIndexes <- which(folds==fold, arr.ind=TRUE)
          
          train_data = filter_features_sets(features_data2[!(features_data2$Uid %in% user_ids[testIndexes]),], !cumul_across, window, slice)
          train_data[,dependent[dep]] = factor(train_data[,dependent[dep]], levels = c(1, 0))
          
          test_data = filter_features_sets(features_data2[(features_data2$Uid %in% user_ids[testIndexes]),], !cumul_across, window, slice)
          test_data[,dependent[dep]] = factor(test_data[,dependent[dep]], levels = c(1, 0))
          
          # Feature selection
          if(feat_selection) {
            removeZeros = apply(train_data, 2,
                                function(x) length(unique(x)) == 1) # remove zeros
            # Don't remove dependent even if only 1 value is available.
            removeZeros[[length(removeZeros)]] = FALSE
            train_data = train_data[, !removeZeros]
            test_data = test_data[, !removeZeros]
            data_cor = cor(train_data[, colnames(train_data) != dependent[dep]])
            data_cor[is.na(data_cor)] = 1
            
            corfeat = findCorrelation(data_cor, cutoff = .85)
            
            if(length(corfeat) > 0) {
              train_data =  train_data[,-corfeat]
              test_data = test_data[,-corfeat]
            }
            fmla <- as.formula(paste0(dependent[dep],"~.", sep=""))
          } else {
            fmla <- as.formula(paste0(dependent[dep],"~", paste0(colnames(train_data[, colnames(train_data) != dependent[dep]]), collapse= "+", sep=""), sep=""))
          }
          
          
          outerloop_train_sets[[fold]] <- train_data
          outerloop_test_sets[[fold]] <- test_data
          
          #==========================================
          
          ctrl <- trainControl(method = "cv", number = 10, allowParallel=T)   # for inner loop
          
          #-------------------------------------------------------------------
          # Grids
          
          rf_Grid   	 <- expand.grid(mtry = c(2,21,40,59,78,98,117,136,155,175))
          
          BstLm_Grid   <- expand.grid(nu=0.1, mstop = seq(5,50,by=5)*10)
          
          xgbTree_Grid <- expand.grid(nrounds = c(50,100), 
                                      max_depth = c(1,2),
                                      eta = c(0.3,0.4), 
                                      gamma = 0,
                                      colsample_bytree = c(0.8,0.6),
                                      min_child_weight = 1,
                                      subsample = c(0.5,1)
                                      )
          
          #-------------------------------------------------------------------------------------------------------------------------
          all_models = c( 
                          "adaboost", "AdaBag", "treebag", "bagFDA", "logicBag", "bagEarth", "gamboost", "glmboost", 
                          "bag", "bartMachine", "bayesglm", "LogitBoost", "J48", "C5.0", "rpartScore", "chaid", "deepboost",
                          "dda", "dwdPoly", "dwdRadial", "RFlda", "fda", "SLAVE", "gpls", "protoclass", "hda", "hdda", "hdrda", 
                          "svmLinearWeights2", "lvq", "lssvmLinear", "lssvmPoly", "lssvmRadial", "lda", "lda2", "stepLDA",
                          "dwdLinear", "svmLinearWeights", "loclda", "LMT", "Mlda", "mda", "manb", "mlpKerasDropoutCost", 
                          "mlpKerasDecayCost", "naive_bayes", "nb", "nbDiscrete", "awnb", "pam", "ORFlog", "ORFpls", "ORFridge",
                          "ORFsvm", "ownn", "polr", "PRIM", "pda", "pda2", "PenalizedLDA", "plr", "multinom", "ordinalNet", "qda", 
                          "stepQDA", "rFerns", "ordinalRF", "rda", "rlda", "regLogistic", "Linda", "rmda", "QdaCov", "rrlda", 
                          "RSimca", "rocc", "rotationForest", "rotationForestCp", "JRip", "PART", "nbSearch", "sda", "CSimca", 
                          "C5.0Rules", "C5.0Tree", "OneR", "sdwd", "sparseLDA", "smda", "slda", "snn", "svmRadialWeights", "tan", 
                          "tanSearch", "awtan", "vbmpRadial", "wsrf", "BstLm", "bstSm", "blackboost", "bstTree", "cforest", "ctree", 
                          "ctree2", "xgbDART", "xgbLinear", "xgbTree", "elm", "gaussprLinear", "gaussprPoly", "gaussprRadial", 
                          "bam", "gam", "gamSpline", "glm", "glmStepAIC", "glmnet", "glmnet_h2o", "gbm_h2o", "kknn", "knn", 
                          "gamLoess", "svmLinear3", "logreg", "avNNet", "monmlp", "mlp", "mlpWeightDecay", "mlpWeightDecayML", "mlpML", 
                          "msaenet", "mlpSGD", "mlpKerasDropout", "mlpKerasDecay", "earth", "gcvEarth", "mxnet", "mxnetAdam", 
                          "nnet", "pcaNNet", "null", "parRF", "partDSA", "kernelpls", "pls", "simpls", "widekernelpls", "plsRglm", 
                          "rbf", "rbfDDA", "ranger", "Rborist", "rf", "extraTrees", "rfRules", "RRF", "RRFglobal", "xyf", "spls", 
                          "dnn", "gbm", "svmBoundrangeString", "svmExpoString", "svmLinear", "svmLinear2", "svmPoly", "svmRadial", 
                          "svmRadialCost", "svmRadialSigma", "svmSpectrumString", "evtree", "nodeHarvest"
                          )
          #-------------------------------------------------------------------------------------------------------------------------
          
          models_names = c("rf", "xgbTree", "BstLm")
          models_grids = list(rf_Grid, xgbTree_Grid, BstLm_Grid)
          
          m = 3

          set.seed(1)
          
          fit <- train(fmla,
                       data = train_data,
                          # method = all_models[],
                       method = models_names[m],
                          # tuneLength = 1,
                       tuneGrid = models_grids[[m]],
                       trControl = ctrl,
                       na.action=na.pass,
                       metric="Accuracy")
          #    print(fit)
          
          #-------------------------------------Save Results----------------------------------------------------
          if(fold==1){
            All_ResultsDF <- (fit$results[,!(names(fit$results) %in% c("Kappa","AccuracySD","KappaSD"))])
            Sum = (fit$results$Accuracy)
            colnames(All_ResultsDF)[which(names(All_ResultsDF) == "Accuracy")] <- paste("Accuracy_", fold, sep = "")
            
          }else{
            Accuracy = fit$results$Accuracy
            All_ResultsDF <- cbind(All_ResultsDF, Accuracy)
            Sum <- Sum + Accuracy
            colnames(All_ResultsDF)[which(names(All_ResultsDF) == "Accuracy")] <- paste("Accuracy_", fold, sep = "")
          }
          
        }
        
        Avg <- Sum / K
        All_ResultsDF <- cbind(All_ResultsDF, Avg)
        All_Results_Trim <- All_ResultsDF[-(ncol(All_ResultsDF)-c(1:K))]
        
        Best <- (All_ResultsDF[All_ResultsDF$Avg==max(All_ResultsDF$Avg),])
        Best <- Best[-(ncol(Best)-c(1:K))]
        # print(Best)  
        
        #-------------------------------------------Generalized Accuracy---------------------------------------------------------
        # outer loop (evaluation):
        
        perform_evaluation = T          # Set to True after the best model has been selected. maually change the model's grid to the optimal for evaluation
                                        # False while searching for the best model (only for efficiency)
        if(perform_evaluation){
          for(fold in 1:K) {
            # print(fold)
            
            # load data for the outer loop evaluation
            training_data <- outerloop_train_sets[[fold]]
            testing_data <- outerloop_test_sets[[fold]]
  
            fmla <- as.formula(paste0(dependent[dep],"~", paste0(colnames(training_data[, colnames(training_data) != dependent[dep]]), collapse= "+", sep=""), sep=""))
            
            ctrl <- trainControl(method = "none", allowParallel=T)    # for outer loop
            
            
            #-------------------------------------------------------------------
            # Best Grids (set after they have been determined then run the outer loop for evaluation)
                
            rf_Grid   	 <- expand.grid(mtry = 21)
                
            BstLm_Grid   <- expand.grid(nu = 0.1, mstop = 200)
                
            xgbTree_Grid <- expand.grid(nrounds = 50, 
                                        max_depth = 1,
                                        eta = 0.3, 
                                        gamma = 0,
                                        colsample_bytree = 0.8,
                                        min_child_weight = 1,
                                        subsample = 1
                                      )
                
            #-------------------------------------------------------------------
            models_names = c("rf", "xgbTree", "BstLm")
            models_grids = list(rf_Grid, xgbTree_Grid, BstLm_Grid)
            #-------------------------------------------------------------------
            m = 3 # choose the intended model 
            
            model_name = models_names[m]
            tuneparam = models_grids[[m]]
  
            set.seed(1)
            models_stats = train_test_model(model_name, tuneparam, models_stats, feature_imp_fname, fold)
            
          }
        }
        print(c("Best Model: ", Best))
        GE_Acc = mean(models_stats[,paste(model_name,"_Accuracy", sep = "")])
        print(c("Generalized Accuracy of the Best Model: ", GE_Acc))
        
        
        #end CrossValidation run
        #---------------------------------------------------------------------------------------------------
        models_stats[is.na(models_stats)] = 0
        
        for (model in models) {
          
          if (model != "MajorityClass"){        # || numRepeat == 1 ){
            
          result = paste(dependent[dep], window, slice, model, numRepeat, sep=",")
          for (metric in metrics) { # for each of the metrics calculate mean and SD and append to create a row for the model
            if (model != "MajorityClass" && metric == "AccuracyHigh") {
              result = weighted_mean_sd(model, "High", models_stats)
            }
            else if (model != "MajorityClass" && metric == "AccuracyLow") {
              result = weighted_mean_sd(model, "Low", models_stats)
            }
            else if (!(model == "MajorityClass" && (metric == "Kappa" || metric == "AccuracyHigh" || metric == "AccuracyLow"))) {
              model_res = models_stats[, paste(model, "_", metric, sep="")] # take the model_metric column from models_stats
              result = paste(result, mean(model_res), sep=",")
              result = paste(result, sd(model_res), sep=",")
              if (metric == "Accuracy" && model != "MajorityClass"){
                print(mean(model_res))
                Sum = Sum + mean(model_res)
              }
              
            } else{
              result = paste(result, 0, sep=",")
              result = paste(result, 0, sep=",")
            }
          }
          sink(out, append = TRUE)
          print(result)
          sink()
        }}
      } #end CrossValidation all repeats
      
    }#end dependant
  }
}

                                        # Create CSV log:
#--------------------------------------------------------------------------------------------------
write.csv(All_ResultsDF,"C:\\Users\\Alire\\Desktop\\All_ResultsDF.csv", row.names = FALSE)
# write.csv(All_Results_Trim,"C:\\Users\\Alire\\Desktop\\All_Results_Trim.csv", row.names = FALSE)
# write.csv(models_stats,"C:\\Users\\Alire\\Desktop\\models_stats.csv", row.names = FALSE)
#--------------------------------------------------------------------------------------------------

stopCluster(cl)

