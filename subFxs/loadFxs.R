
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
    junk = read.csv(sprintf("%s/%s", dataDir, fileName), col.names = colNames, header = T)
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


