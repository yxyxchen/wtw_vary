
# loadData.R
# varying-magnitude WTW

# each subject has:
#  header file (*hdr.txt)
#  main data file (ending both in .mat and .txt)
#  empty *keytimes.txt file, a vestige of the effort key-pressing version. 
loadAllData = function() {
  library("gtools")
  # load hdrData fileNames
  dataDir = "data"
  fileNames = list.files(path= dataDir, pattern=('*_hdr.txt'))
  fileNames = mixedsort(sort(fileNames))
  nFile = length(fileNames)
  if(any(duplicated(fileNames))){
    print("duplicated hdrData files!")
    break
  }else{
    sprintf("load %d hdrData files", nFile) 
  }
  # load hdrData
  ID = vector(length = nFile)
  cbal = vector(length = nFile)
  for(i in 1 : nFile){
    fileName = fileNames[i]
    id = substr(fileName, 12,14)
    junk = read.csv(sprintf("%s/%s", dataDir, fileName))
    ID[i] = id
    cbal[i] = as.double(substr(junk[[1]][1], 7, 7))
  }
  hdrData = data.frame(ID = ID, cbal = cbal, stringsAsFactors = F)
  
  # define column names 
  colNames = c('blockNum', 'trialNum', 'trialStartTime', 'nKeyPresses', 'scheduledWait',
    'rewardTime', 'timeWaited', 'sellTime', 'trialEarnings','totalEarnings')
  
  # get trialData fileNames
  trialFileNames = list.files(path= dataDir, pattern=('wtw-work-7_[0-9]{3}_[0-9].txt'))
  trialFileNames = mixedsort(sort(trialFileNames))
  nFile = length(trialFileNames)
  if(any(duplicated(trialFileNames))){
    print("duplicated trialData files!")
    break
  }else{
    sprintf("load %d trialData files", nFile) 
  }
  
  # load trialData
  trialData = list()
  for(i in 1 : nFile){
    fileName = trialFileNames[i]
    id = substr(fileName, 12,14)
    junk = read.csv(sprintf("%s/%s", dataDir, fileName), col.names = colNames)
    if(hdrData$cbal[which(hdrData$ID == id)] == 1){
      condition = ifelse(junk$blockNum %% 2 == 1, "Rising", "Falling")
    }else{
      condition = ifelse(junk$blockNum %% 2 == 1, "Falling", "Rising")
    }
    junk$condition = condition
    trialData[[id]] = junk
  }
  
  outputData = list(hdrData=hdrData, trialData=trialData)
  return(outputData)
  
} 


loadExpPara = function(paras, dirName){
  # number of paras 
  nE = length(paras) + 1
  # number of files
  fileNames = list.files(path= dirName, pattern=("*_summary.txt"))
  library("gtools")
  fileNames = mixedsort(sort(fileNames))
  n = length(fileNames) 
  sprintf("load %d files", n)
  
  # initialize the outout variable 
  expPara = matrix(NA, n, nE * 4)
  idList = vector(length = n)
  # loop over files
  for(i in 1 : n){
    fileName = fileNames[[i]]
    address = sprintf("%s/%s", dirName, fileName)
    junk = read.csv(address, header = F)
    idIndexs = str_locate(fileName, "s[0-9]+")
    idList[i] = as.double(substr(fileName, idIndexs[1]+1, idIndexs[2]))
    # delete the lp__ in the old version
    if(nrow(junk) > nE){
      junk = junk[1:nE,]
    }
    expPara[i, 1:nE] = junk[,1]
    expPara[i, (nE + 1) : (2 * nE)] = junk[,2]
    expPara[i, (2*nE + 1) : (3 * nE)] = junk[,9]
    expPara[i, (3 * nE + 1) : (4 * nE)] = junk[,10]
  }
  # transfer expPara to data.frame
  expPara = data.frame(expPara)
  junk = c(paras, "LL_all")
  colnames(expPara) = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"))
  expPara$id = idList
  return(expPara)
}

getMode = function(x){
  ob = density(x)
  value = ob$x
  dense = ob$y
  return(value[which.max(dense)])
}

loadExpParaExtra = function(modelName, paras){
  nE = length(paras) + 2
  load("genData/expDataAnalysis/blockData.RData")
  idList = unique(blockData$id) 
  n = length(idList)
  nE = length(paras) + 1
  expParaMode = matrix(NA, n, nE)
  expParaMedian = matrix(NA, n, nE)
  for(i in 1 : n){
    ID = idList[i]
    fileName = sprintf("genData/expModelFitting/%s/s%d.txt", modelName, ID)
    junk = read.csv(fileName, header = F)
    
    expParaMode[i, ] = apply(junk[,1:nE], MARGIN = 2, getMode)
    expParaMedian[i, ] = apply(junk[,1:nE], MARGIN = 2, median)
  }
  expParaMode = data.frame(expParaMode)
  expParaMedian = data.frame(expParaMedian)
  
  junk = c(paras, "log_lik")
  colnames(expParaMode) =  paste0(junk, "Mode")
  expParaMode$id = idList  # needed for left_join
  colnames(expParaMedian) =  paste0(junk, "Median")
  expParaMedian$id = idList  # needed for left_join
  
  outputs = list("expParaMedian" = expParaMedian, "expParaMode" = expParaMode)
  return(outputs)
}

# for each sub
loadExpParaMedian = function(modelName, paras, id){
  nE = length(paras) + 2
  load("genData/expDataAnalysis/blockData.RData")
  nE = length(paras) + 1
  expParaMedian = vector(length = nE)
    fileName = sprintf("genData/expModelFitting/%s/s%d.txt", modelName, id)
    junk = read.csv(fileName, header = F)
    expParaMedian  = apply(junk[,1:nE], MARGIN = 2, median)
  return(expParaMedian)
}

loadSimPara = function(paras, dirName){
  # number of paras 
  nE = length(paras) + 1
  # number of files
  fileNames = list.files(path= dirName, pattern=("*_summary.txt"))
  library("gtools")
  fileNames = mixedsort(sort(fileNames))
  n = length(fileNames) 
  sprintf("load %d files", n)
  
  # extract seqIdxs and cbIdxs from filenames
  seqIdxs = as.double(sapply(1:n, function(i){
    tempt = regexpr("_r[0-9]*", fileNames[i]) 
    match.length = attr(tempt, "match.length")
    start = tempt[[1]] + 2
    ending = tempt[[1]] + match.length - 1
    substring(fileNames[i], start, ending)
  })) + 1
  nSeq = length(unique(seqIdxs))
  
  cbIdxs =  as.double(sapply(1:n, function(i){
    tempt = regexpr("_s[0-9]*", fileNames[i]) 
    match.length = attr(tempt, "match.length")
    start = tempt[[1]] + 2
    ending = tempt[[1]] + match.length - 1
    substring(fileNames[i], start, ending)
  }))
  nComb = length(unique(cbIdxs))
  
  # initialize the outout variable 
  junk = c(paras, "LL_all")
  varNames = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"))
  simParaHP = array(dim = c(nE * 4, nSeq, nComb),
                    dimnames = list(varNames, 1 : nSeq, 1 : nComb))
  simParaLP = array(dim = c(nE * 4, nSeq, nComb), dimnames = list(varNames, 1 : nSeq, 1 : nComb))
  # loop over files
  for(i in 1 : n){
    # read the file
    fileName = fileNames[[i]]
    address = sprintf("%s/%s", dirName, fileName)
    junk = read.csv(address, header = F)
    
    # determine the condition 
    if(grepl("HP", fileName)){
      simParaHP[, seqIdxs[i], cbIdxs[i]] = unlist(junk[,c(1,2,9,10)])
    }else{
      simParaLP[, seqIdxs[i], cbIdxs[i]] = unlist(junk[,c(1,2,9,10)])
    }
  }
  simPara = list(HP = simParaHP, LP =simParaLP)
  # transfer expPara to data.frame
  return(simPara)
}

