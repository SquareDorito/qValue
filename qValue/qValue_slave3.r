parseTextFile = function(counters_FilePath){
  link_table = read.table(counters_FilePath,sep='\t',header=FALSE)
  link_df=data.frame(link_table)
  rownames(link_df)<-link_df[,2]
  write.table(link_df,file="link_df.txt",sep="\t",quote=FALSE,row.names=FALSE)
  return(link_df)
}

parseFileHelper = function(filePath, waitPath, replicates, timepts, label_free, unpaired, silac_num, silac_den, min_rep, bh, minPA, maxPA, is_min, log_LF, log_silac, is_error) {
  ## Assert that the XML filePath exists, exit if it doesn't exist.
  print(paste("Parsing ", filePath, sep=""))
  
  parsedXML = parseXML(filePath, timepts, replicates, minPA, maxPA, is_min) 
  ## Creates a copy of the parsed XML matrix to check errors.
  ## write.table(parsedXML, "result.csv", row.names=FALSE, col.names=FALSE, sep=",")
  
  
  ## Turn parsedXML matrix into data frame.
  dim_pXML = dim(parsedXML)
  pXML_pep = parsedXML[1,]
  pXML_data = parsedXML[2:dim_pXML[1],]
  pXML_df = data.frame(pXML_data)
  colnames(pXML_df) = pXML_pep
  p_counters<<-pXML_pep
  rNameCounter = 1
  for (rp in 1:replicates){
    for (ch in 1:3){
      for (tp in 1:timepts){
        tempRN = paste("ch", toString(ch), " rep", toString(rp), " timepoint", toString(tp), sep = "")
        ##print(paste0(tempRN))
        ##tempRN = paste(paste(paste("ch", toString(ch), sep = ""), paste(" rep", toString(rp), sep = "")),paste(" timepoint", toString(tp),sep = ""))
        rownames(pXML_df)[rNameCounter] <- tempRN
        rNameCounter = rNameCounter + 1
      }
    }
  }
  
  
  if(identical("C:\\Users\\knoh1\\Documents\\qValue_lf_test\\phospho.xml",filePath)){
    write.table(pXML_df,file="phospho_df.txt",sep=",",col.names=colnames(pXML_df),row.names=rownames(pXML_df),quote=FALSE)
  }else{
    write.table(pXML_df,file="unphospho_df.txt",sep=",",col.names=colnames(pXML_df),row.names=rownames(pXML_df),quote=FALSE)
  }
  return(pXML_df)
}

parseXML = function(path, timepts, replicates, minPA, maxPA, is_min) {
  ## These variables will be used to create the dimensions of resultmtx.
  print(timepts)
  print(replicates)
  norows = (timepts*replicates*3)+1
  nocols = 0
  ## These counters will help the XML parsing keeps track of where in the matrix to place and store variable.
  rowCounter = 1
  colCounter = 1
  dataCounter = 1
  ## This matrix will be returned as the result of parsing the XML file.
  resultmtx = matrix()
  ## These booleans keep track of COL and DATA elements.
  colBool = FALSE
  datBool = FALSE
  ## NotANumber will be put in place of empty elements in the matrix.
  nan = (1/0)-(1/0)
  xmlEventParse(
    path, 
    handlers = list(
      startElement = function(name, attr) {
        if (name == "FIELD") {
          ## ignore tag
        }
        if (name == "RESULTSET") {
          ## RESULTSET tag will come with a FOUND="##" <RESULTSET FOUND="99">.
          ## That number is the number of peptides, we want the number of columns to be 1 + that number.
          nocols <<- attr[1]
          nocols <<- as.numeric(nocols) + 1
          ## Initialize the dimensions of resultmtx.
          resultmtx <<- matrix(nrow=norows,ncol=nocols)
        }
        if (name == "ROW") {
          ## Set the rowCounter back to 1. Increment the colCounter.
          rowCounter <<- 1
          colCounter <<- colCounter + 1
        }
        if (name == "COL") {
          colBool <<- TRUE
        }
        if (name == "DATA") {
          datBool <<- TRUE
        }
      }, text = function(text) {
        ## <DATA> elements are followed by text holding the value of the data this needs to be compared to the threshold and stored in the resultmtx.
        if (datBool) {
          ## Note that the Peptide Counter is always the first <DATA> element after a new <ROW> element and it's important to not threshold these counter numbers.
          if(rowCounter == 1){
            resultmtx[rowCounter,colCounter] <<- as.numeric(text)
          } else {
            ## Compare the text to the necessary peak area threshold.
            if (is_min){
              if (as.numeric(text) > minPA) {
                resultmtx[rowCounter,colCounter] <<- as.numeric(text)
              }
            } else {
              if (as.numeric(text) < maxPA) {
                resultmtx[rowCounter,colCounter] <<- as.numeric(text)
              }
            }
          }
        }
      }, endElement = function(name, uri) {
        if (name == "COL") {
          ## Set the dataCounter back to 1.
          dataCounter <<- 1
          colBool <<- FALSE
        }
        if (name == "DATA") {
          ## Check the dataCounter to know when to stop reading input.
          ## The input is hard coded to contain 16 timepoints but in practice we often use less.
          ## Any timepoints beyond the number used are set as NA and we want to ignore these values.
          if (dataCounter <= timepts) {
            ## Any thresholded or empty elements in the matrix need to be set as NaN
            if (is.na(resultmtx[rowCounter,colCounter])) {
              resultmtx[rowCounter,colCounter] <<- nan
            }
            rowCounter <<- rowCounter + 1
          }
          dataCounter <<- dataCounter + 1
          datBool <<- FALSE
        }
      }
    ))
  return(resultmtx)
}

create_lf_df = function(p_df,np_df,link_df,replicates,timepts,min_rep,label_free,unpaired,silac_num,silac_den,bh,minPA,maxPA,is_min,log_LF,log_silac,is_error){
  #delete then recreate local_qValue folder on C drive to ensure the folder is empty
  setwd("C:\\")
  fn1<-'local_qValue'
  if(file.exists(fn1)){
    unlink(fn1,recursive=TRUE)
  }
  dir.create(fn1)
  setwd("C:\\local_qValue")
  fn2<-'labelfree_df'
  if(file.exists(fn2)){
    unlink(fn2,recursive=TRUE)
  }
  dir.create(fn2)
  
  p_list<-link_df[,2]
  p_num=length(p_list)
  
  if(testing==1){
    #print out the phospho counters for testing
    p_counters_out<-data.frame(p_list)
    write.table(p_counters_out,file="p_counters_out.txt",sep="",quote=FALSE,row.names=FALSE)
  }
  
  #will store the dataframes in a global list to be used in ANOVA later
  lf_df_vector <<- vector("list",length=timepts*p_num)
  
  empty_numeric<-numeric()
  empty_char<-character()
  
  for(i in 1:p_num){
    for (tp in 1:timepts){
      repcount_a=0
      repcount_b=0
      repcount_c=0
      repcount_d=0
      id<-character()
      np<-numeric()
      ko<-numeric()
      peakArea<-numeric()
      for (rp in 1:replicates){
        counter=toString(p_list[i])
        if(!counter %in% colnames(p_df)){
          next
        }
        #the rownames associated with the silac numerator and den
        tempRN1 = paste("ch", toString(silac_num), " rep", toString(rp), " timepoint", toString(tp), sep = "")
        if(is.na(p_df[tempRN1,counter])){
          next
        }
        if(testing==1){
          cat(paste("p_counter: ",counter," peakArea_num: ",p_df[tempRN1,counter],"\n"))
        }
        repcount_a=repcount_a+1
        tempID1 = paste("pep",counter,"_rep", toString(rp), sep = "")
        ko<-c(ko,0)
        id<-c(id,tempID1)
        np<-c(np,0)
        peakArea<-c(peakArea,p_df[tempRN1,counter])
      }
      for (rp in 1:replicates){
        counter=toString(p_list[i])
        #the rownames associated with the silac numerator and den
        tempRN1 = paste("ch", toString(silac_num), " rep", toString(rp), " timepoint", toString(tp), sep = "")
        tempID1 = paste("pep",counter,"_rep", toString(rp), sep = "")
        rawNP<-as.character(link_df[counter,3])
        #list of np counters associated with the phospho counter
        np_list<-as.list(strsplit(rawNP,",")[[1]])
        if(!length(np_list)){
          next
        }
        repcount_b=repcount_b+1
        for(j in 1:length(np_list)){
          np_counter<-toString(np_list[[j]])
          if(testing==1){
            cat(paste("np_counter: ",np_counter," peakArea_num: ",np_df[tempRN1,np_counter],"\n"))
          }
          if(!np_counter %in% colnames(np_df)){
            next
          }
          if(is.na(np_df[tempRN1,np_counter])){
            next
          }
          ko<-c(ko,0)
          id<-c(id,tempID1)
          np<-c(np,1)
          peakArea<-c(peakArea,np_df[tempRN1,np_counter])
        }
      }
      for (rp in 1:replicates){
        counter=toString(p_list[i])
        if(!counter %in% colnames(p_df)){
          next
        }
        #the rownames associated with the silac numerator and den
        tempRN2 = paste("ch", toString(silac_den), " rep", toString(rp), " timepoint", toString(tp), sep = "")
        if(is.na(p_df[tempRN2,counter])){
          next
        }
        if(testing==1){
          cat(paste("p_counter: ",counter," peakArea_num: ",p_df[tempRN2,counter],"\n"))
        }
        repcount_c=repcount_c+1
        tempID2 = paste("pep",counter,"_rep", toString(rp), sep = "")
        ko<-c(ko,1)
        id<-c(id,tempID2)
        np<-c(np,0)
        peakArea<-c(peakArea,p_df[tempRN2,counter])
      }
      for (rp in 1:replicates){
        counter=toString(p_list[i])
        if(!counter %in% colnames(p_df)){
          next
        }
        #the rownames associated with the silac numerator and den
        tempRN2 = paste("ch", toString(silac_den), " rep", toString(rp), " timepoint", toString(tp), sep = "")
        tempID2 = paste("pep",counter,"_rep", toString(rp), sep = "")
        rawNP<-as.character(link_df[counter,3])
        #list of np counters associated with the phospho counter
        np_list<-as.list(strsplit(rawNP,",")[[1]])
        if(!length(np_list)){
          next
        }
        repcount_d=repcount_d+1
        for(j in 1:length(np_list)){
          np_counter<-toString(np_list[[j]])
          if(testing==1){
            cat(paste("np_counter: ",np_counter," peakArea_num: ",np_df[tempRN2,np_counter],"\n"))
          }
          if(!np_counter %in% colnames(np_df)){
            next
          }
          if(is.na(np_df[tempRN2,np_counter])){
            next
          }
          ko<-c(ko,1)
          id<-c(id,tempID2)
          np<-c(np,1)
          peakArea<-c(peakArea,np_df[tempRN2,np_counter])
        }
      }
      #storing the generated dataframe into global list
      temp_index=(timepts*(i-1))+tp
      if(repcount_a>=min_rep && repcount_b>=min_rep && repcount_c>=min_rep && repcount_d>=min_rep){
        lf_df_vector[[temp_index]]<<-data.frame(id,ko,np,peakArea,row.names=NULL)
      }else{
        lf_df_vector[[temp_index]]<<-data.frame(empty_char,empty_numeric,empty_numeric,empty_numeric,row.names=NULL)
      }
      
      if(testing==1){
        #print the dataframe to folder for testing purposes
        outfile=paste(fn2,"/sample",toString(p_list[i]),"_timepoint",toString(tp),".txt",sep="")
        write.table(lf_df_vector[[temp_index]],file=outfile,sep="\t",row.names=FALSE,quote=FALSE)
      }
      rm(id)
      rm(np)
      rm(ko)
      rm(peakArea)
    }
  }
}

create_silac_df = function(p_df,np_df,link_df,replicates,timepts,min_rep,label_free,unpaired,silac_num,silac_den,bh,minPA,maxPA,is_min,log_LF,log_silac,is_error,testing){
  
  setwd("C:\\")
  fn1<-'local_qValue'
  if(file.exists(fn1)){
    unlink(fn1,recursive=TRUE)
  }
  dir.create(fn1)
  setwd("C:\\local_qValue")
  fn2<-'silac_df'
  if(file.exists(fn2)){
    unlink(fn2,recursive=TRUE)
  }
  dir.create(fn2)
  
  empty_numeric<-numeric()
  empty_char<-character()
  
  p_list<-link_df[,2]
  p_num=length(p_list)
  
  if(testing==1){
    p_counters_out<-data.frame(p_list)
    write.table(p_counters_out,file="p_counters_out.txt",sep="",quote=FALSE,row.names=FALSE)
  }
  
  silac_df_vector <<- vector("list",length=timepts*p_num)
  
  for(i in 1:p_num){
    for (tp in 1:timepts){
      id<-character()
      np<-numeric()
      ratio<-numeric()
      repcount_p=0
      repcount_np=0
      for (rp in 1:replicates){
        counter=toString(p_list[i])
        num_rn = paste("ch", toString(silac_num), " rep", toString(rp), " timepoint", toString(tp), sep = "")
        den_rn = paste("ch", toString(silac_den), " rep", toString(rp), " timepoint", toString(tp), sep = "")
        if(testing==1){
          cat(paste("p_counter: ",counter," ",num_rn," ",p_df[num_rn,counter],"\n",sep=""))
          cat(paste("p_counter: ",counter," ",den_rn," ",p_df[den_rn,counter],"\n",sep=""))
        }
        if(is.na(p_df[num_rn,counter]) || is.na(p_df[den_rn,counter])){
          next
        }
        repcount_p=repcount_p+1
        tempID = paste("pep",counter,"_rep", toString(rp), sep = "")
        tempRatio=(p_df[num_rn,counter])/(p_df[den_rn,counter])
        if(testing==1){
          cat(paste("Ratio: ",tempRatio,"\n",sep=""))
        }
        id<-c(id,tempID)
        np<-c(np,0)
        ratio<-c(ratio,tempRatio)
        rawNP<-as.character(link_df[counter,3])
        np_list<-as.list(strsplit(rawNP,",")[[1]])
        if(!length(np_list)){
          next
        }
        for(j in 1:length(np_list)){
          np_counter<-np_list[[j]]
          if(testing==1){
            cat(paste("np_counter: ",np_counter," peakArea_num: ",np_df[num_rn,np_counter],"\n"))
            cat(paste("np_counter: ",np_counter," peakArea_den: ",np_df[den_rn,np_counter],"\n"))
          }
          if(!np_list[[j]] %in% colnames(np_df)){
            next
          }
          if(is.na(np_df[num_rn,np_list[[j]]])){
            next
          }
          if(is.na(np_df[den_rn,np_list[[j]]])){
            next
          }
          repcount_np=repcount_np+1
          id<-c(id,tempID)
          np<-c(np,1)
          np_ratio=(np_df[num_rn,np_list[[j]]])/(np_df[den_rn,np_list[[j]]])
          ratio<-c(ratio,np_ratio)
        }
      }
      temp_index=(timepts*(i-1))+tp
      if(repcount_p>=min_rep && repcount_np>=min_rep){
        silac_df_vector[[temp_index]]<<-data.frame(id,np,ratio)
      } else{
        silac_df_vector[[temp_index]]<<-data.frame(empty_char,empty_numeric,empty_numeric)
      }
      if(testing==1){
        outfile=paste(fn2,"/sample",p_list[i],"_timepoint",toString(tp),".txt",sep="")
        write.table(silac_df_vector[[temp_index]],file=outfile,sep="\t",row.names=FALSE,quote=FALSE)
      }
      rm(id)
      rm(np)
      rm(ratio)
    }
  }
}

calc_qVals=function(p_df,np_df,link_df,timepts,unpaired,bh,testing){
  library(qvalue)
  setwd("C:\\local_qValue")
  stats_df<-data.frame()
  fileNames<-character()
  pVals<-numeric()
  qVals<-numeric()
  p_list<-link_df[,2]
  no_qVals_matrix<-matrix(0,nrow=length(p_list),ncol=timepts)
  no_qVals<-numeric()
  
  for(i in 1:length(p_list)){
    for(j in 1:timepts){
      if(!unpaired){
        #silac stuff here
        tempFileName=paste("silac_df/","sample",toString(p_list[i]),"_timepoint",j,".txt",sep="")
        fileNames<-c(fileNames,tempFileName)
        if(!file.exists(tempFileName)){
          no_qVals_matrix[i,j]=1
          next
        }
        temp_index=(timepts*(i-1))+j
        temp_df=silac_df_vector[[temp_index]]
        dat<-temp_df
        
        if(is.data.frame(dat) && nrow(dat)==0){
          no_qVals_matrix[i,j]=1
          next
        }
        model = lm(log(ratio) ~ np, data = dat, na.action = na.omit)
        a=anova(model)
        if(is.na(a[['Pr(>F)']][1])){
          no_qVals_matrix[i,j]=1
        }else{
          pVals<-c(pVals,a[['Pr(>F)']][1])
        }
      }else{
        #label free stuff
        tempFileName=paste("labelfree_df/","sample",colnames(p_df)[i],"_timepoint",j,".txt",sep="")
        fileNames<-c(fileNames,tempFileName)
        
        if(!file.exists(tempFileName)){
          no_qVals_matrix[i,j]=1
          next
        }
        temp_index=(timepts*(i-1))+j
        temp_df=lf_df_vector[[temp_index]]
        dat<-temp_df
        
        if(is.data.frame(dat) && nrow(dat)==0){
          no_qVals_matrix[i,j]=1
          next
        }
        
        model1 = lm(log(peakArea) ~ ko + np, data = dat, na.action = na.omit)
        model2 = lm(log(peakArea) ~ ko*np, data = dat, na.action = na.omit)
        a = anova(model1, model2)
        if(is.na(a[['Pr(>F)']][2])){
          no_qVals_matrix[i,j]=1
        }else{
          pVals<-c(pVals,a[['Pr(>F)']][2])
        }
      }
    }
  }
  
  p_vals_out<-data.frame(pVals)
  write.table(p_vals_out,file="pvals.txt",sep="\t",quote=FALSE,row.names=FALSE)
  
  if(bh==1){
    print(paste("Controlling FDR with the Benjamin-Hochberg method.\n",sep=""))
    qobj<-qvalue(p=pVals,lambda=0)
  }else{
    print(paste("Controlling FDR with the Storey method.\n",sep=""))
    qobj<-qvalue(p=pVals)
  }
  qvalues <- qobj$qvalues
  if(length(qvalues)==0){
    qobj <- qvalue(pvals, pi0.method="bootstrap")
    qvalues <- qobj$qvalues
  }
  upper=0.95
  while(length(qvalues)==0){
    qobj <- qvalue(pvals, lambda=seq(0, " + upper + ", .01))
    qvalues <- qobj$qvalue
    upper=upper-0.01
    if(upper<0.01){
      break
    }
  }
  
  lower=0
  while(length(qvalues)==0){
    qobj <- qvalue(pvals, lambda=seq(" + lower + ", .99, .01))
    qvalues<-qobj$qvalues
    lower=lower+0.01
    if(lower>0.98){
      break
    }
  }
  
  if(length(qvalues)==0){
    qobj<-qvalue(pvals,lambda=0)
    qvalues<-qobj$qvalues
  }
  
  qVals<-c(qVals,qvalues)
  output_qVals<-matrix(,nrow=length(p_list),ncol=2*timepts)
  norm_counter_list<-link_df[,1]
  
  qVal_counter=1
  
  for(i in 1:length(p_list)){
    for(j in 1:timepts){
      if(no_qVals_matrix[i,j]==0){
        output_qVals[i,j]=qVals[qVal_counter]
        output_qVals[i,timepts+j]=pVals[qVal_counter]
        qVal_counter=qVal_counter+1
      }
    }
  }
  
  write.table(output_qVals,file="output.txt",sep="\t",row.names=norm_counter_list,col.names=FALSE,quote=FALSE)
  
  phist<-"phist.png"
  png(filename=phist)
  hist(pVals,breaks=100)
  dev.off()
  qplot<-"qplot.png"
  png(filename=qplot)
  qplot=plot(qobj)
  dev.off()
  data<-cbind(pVals,qVals)
  if(testing==1){
    return(qVals)
  }
}

library(XML)
args = commandArgs(trailingOnly =TRUE)
## The following are arguments thatm ust be passed
p_filePath = "NA"
np_filePath = "NA"
waitPath = "NA"
counters_filePath="NA"
replicates = -1
timepts = -1
min_rep = 3
## The following are optional arguments that may be passed
label_free = 2
unpaired = 0
silac_num = 2
silac_den = 3
bh = 1
minPA = -Inf
maxPA = Inf
is_min = TRUE
log_LF = FALSE
log_silac = FALSE
is_error = FALSE
testing=0

if (length(args) < 7) {
  stop("Invalid Input: Too Few Arguments.", call.=FALSE)
} else if (length(args) > 18) {
  stop("Invalid Input: Too Many Arguments.", call.=FALSE)
} else {
  ## The first two arguments when called must be a path to the XML file and a path to the wait file.
  p_filePath = args[1]
  np_filePath = args[2]
  counters_filePath = args[3]
  waitPath = args[4]
  ## This for loop goes through the rest of the arguments and splits each element into a 1x2 matrix where the elements are
  ## the sides of the equals sign (ex. arg_mtx = [timepoints, 12]).
  ## Note that for arguments without an equals sign each element in the matrix will be the same (ex. arg_mtx = [min, min]).
  for (i in 5:length(args)){
    arg_mtx=matrix(strsplit(args[i], "[=]")[[1]], nrow=1, ncol=2)
    if (arg_mtx[1,1] == "replicates"){
      replicates = as.numeric(arg_mtx[1,2])
    } else if (arg_mtx[1,1] == "timepoints"){
      timepts = as.numeric(arg_mtx[1,2])
    } else if (arg_mtx[1,1] == "label_free"){
      label_free = as.numeric(arg_mtx[1,2])
    } else if (arg_mtx[1,1] == "unpaired"){
      unpaired = as.numeric(arg_mtx[1,2])
    } else if (arg_mtx[1,1] == "silac_numerator"){
      silac_num = as.numeric(arg_mtx[1,2])
    } else if (arg_mtx[1,1] == "silac_denominator"){
      silac_den = as.numeric(arg_mtx[1,2])
    } else if (arg_mtx[1,1] == "min_replicates"){
      min_rep = as.numeric(arg_mtx[1,2])
    } else if (arg_mtx[1,1] == "bh"){
      bh = as.numeric(arg_mtx[1,2])
    } else if (arg_mtx[1,1] == "min_peak_area"){
      ## Note that is_min is checked here meaning that min/max should be specified as an argument before min_peak_area.
      if (is_min) {
        minPA = as.numeric(arg_mtx[1,2])
      } else {
        maxPA = as.numeric(arg_mtx[1,2])
      }
    } else if (arg_mtx[1,1] == "min"){
      is_min = TRUE
    } else if (arg_mtx[1,1] == "max"){
      is_min = FALSE
    } else if (arg_mtx[1,1] == "logLF"){
      log_LF = TRUE
    } else if (arg_mtx[1,1] == "logSILAC"){
      log_silac = TRUE
    } else if (arg_mtx[1,1] == "error"){
      is_error = TRUE
    } else if (arg_mtx[1,1] == "testing"){
      testing = as.numeric(arg_mtx[1,2])
    } else {
      stop("Invalid Input: Unknown Argument.", call.=FALSE)
    }
  }
  p_df<-parseFileHelper(p_filePath, waitPath, replicates, timepts, label_free, unpaired, silac_num, silac_den, min_rep, bh, minPA, maxPA, is_min, log_LF, log_silac, is_error, testing)
  #setwd('C:\\Users\\knoh1\\Documents\\qValue_product')
  #np_df = read.table('unphospho_df.txt',sep=',',header=TRUE,check.names=FALSE)
  np_df<-parseFileHelper(np_filePath, waitPath, replicates, timepts, label_free, unpaired, silac_num, silac_den, min_rep, bh, minPA, maxPA, is_min, log_LF, log_silac, is_error, testing)
  link_df<-parseTextFile(counters_filePath)
  if(unpaired==0){
    create_silac_df(p_df,np_df,link_df,replicates,timepts,min_rep,label_free,unpaired,silac_num,silac_den,bh,minPA,maxPA,is_min,log_LF,log_silac,is_error, testing)
  } else{
    create_lf_df(p_df,np_df,link_df,replicates,timepts,min_rep,label_free,unpaired,silac_num,silac_den,bh,minPA,maxPA,is_min,log_LF,log_silac,is_error, testing)
  }
  calc_qVals(p_df,np_df,link_df,timepts,unpaired,bh,testing)
}

cat("Calculations completed.")