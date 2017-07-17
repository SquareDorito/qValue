parseTextFile = function(){
  setwd('C:/Users/knoh1/Documents/qValue_product')
  link_table = read.table('p_to_np.txt',sep='\t',header=TRUE)
  link_df=data.frame(link_table)
  rownames(link_df)<-link_df[,2]
  write.table(link_df,file="link_df.txt",sep="\t",quote=FALSE)
  return(link_df)
}

parseFileHelper = function(filePath, waitPath, replicates, timepts, label_free, unpaired, silac_num, silac_den, min_rep, bh, minPA, maxPA, is_min, log_LF, log_silac, is_error) {
    ## Assert that the XML filePath exists, exit if it doesn't exist.
	print(paste("Parsing ", filePath, sep=""))
    if (!file.exists(filePath)) waitForInput(paste(filePath, " does not exist. Please check input folder.", sep=""))

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
	#cat(pXML_pep,file="colnames.txt",sep="\n")
	rNameCounter = 1
	for (ch in 1:3){
		for (rp in 1:replicates){
			for (tp in 1:timepts){
				tempRN = paste("ch", toString(ch), " rep", toString(rp), " timepoint", toString(tp), sep = "")
				##print(paste0(tempRN))
				##tempRN = paste(paste(paste("ch", toString(ch), sep = ""), paste(" rep", toString(rp), sep = "")),paste(" timepoint", toString(tp),sep = ""))
				rownames(pXML_df)[rNameCounter] <- tempRN
				rNameCounter = rNameCounter + 1
			}
		}
	}
	
	if(identical("C:\\Users\\knoh1\\Documents\\qValue_product\\phospho_data.xml",filePath)){
	  write.table(pXML_df,file="phospho_df.txt",sep=",",col.names=colnames(pXML_df),row.names=rownames(pXML_df),quote=FALSE)
	}else{
	  write.table(pXML_df,file="unphospho_df.txt",sep=",",col.names=colnames(pXML_df),row.names=rownames(pXML_df),quote=FALSE)
	}
	return(pXML_df)
}

parseXML = function(path, timepts, replicates, minPA, maxPA, is_min) {
    ## These variables will be used to create the dimensions of resultmtx.
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
  
  setwd('C:/Users/knoh1/Documents/qValue_product')
  fn<-'labelfree_df'
  if(file.exists(fn)){
    unlink(fn,recursive=TRUE)
  }
  dir.create(fn)
  
  p_list<-link_df[,2]
  p_num=length(p_list)
  
  for(i in 1:p_num){
    for (tp in 1:timepts){
      answer <- data.frame()
      id<-character()
      np<-numeric()
      channel<-numeric()
      peakArea<-numeric()
      for (ch in 1:3){
        for (rp in 1:replicates){
          counter=toString(p_list[i])
          tempRN = paste("ch", toString(ch), " rep", toString(rp), " timepoint", toString(tp), sep = "")
          if(is.na(p_df[tempRN,counter])){
            next
          }
          tempID = paste("pep",counter,"_rep", toString(rp), sep = "")
          channel<-c(channel,ch)
          id<-c(id,tempID)
          np<-c(np,0)
          peakArea<-c(peakArea,p_df[tempRN,counter])
          rawNP<-as.character(link_df[counter,3])
          np_list<-as.list(strsplit(rawNP,",")[[1]])
          if(!length(np_list)){
            next
          }
          for(j in 1:length(np_list)){
            if(is.na(np_df[tempRN,np_list[[j]]])){
              next
            }
            channel<-c(channel,ch)
            id<-c(id,tempID)
            np<-c(np,1)
            peakArea<-c(peakArea,np_df[tempRN,np_list[[j]]])
          }
        }
      }
      if(length(peakArea)!=0){
        answer <- data.frame(id,channel,np,peakArea)
        outfile=paste(fn,"/sample",toString(p_list[i]),"_timepoint",toString(tp),".txt",sep="")
        write.table(answer,file=outfile,sep="\t",row.names=FALSE,quote=FALSE)
      }
      rm(answer)
      rm(id)
      rm(np)
      rm(channel)
      rm(peakArea)
    }
  }
}

create_silac_df = function(p_df,np_df,link_df,replicates,timepts,min_rep,label_free,unpaired,silac_num,silac_den,bh,minPA,maxPA,is_min,log_LF,log_silac,is_error){
  
  setwd('C:/Users/knoh1/Documents/qValue_product')
  fn<-'silac_df'
  if(file.exists(fn)){
    unlink(fn,recursive=TRUE)
  }
  dir.create(fn)
  
  p_list<-link_df[,2]
  p_num=length(p_list)
  
  for(i in 1:p_num){
    for (tp in 1:timepts){
      answer <- data.frame()
      id<-character()
      np<-numeric()
      ratio<-numeric()
      for (rp in 1:replicates){
        counter=toString(p_list[i])
        num_rn = paste("ch", toString(silac_num), " rep", toString(rp), " timepoint", toString(tp), sep = "")
        den_rn = paste("ch", toString(silac_den), " rep", toString(rp), " timepoint", toString(tp), sep = "")
        if(is.na(p_df[num_rn,counter]) || is.na(p_df[den_rn,counter])){
          next
        }
        tempID = paste("pep",counter,"_rep", toString(rp), sep = "")
        tempRatio=(p_df[num_rn,counter])/(p_df[den_rn,counter])
        #cat(paste(tempRatio,"\n",sep=""))
        id<-c(id,tempID)
        np<-c(np,0)
        ratio<-c(ratio,tempRatio)
        if(!any(row.names(link_df) == counter)){
          next
        }
        rawNP<-as.character(link_df[counter,3])
        np_list<-as.list(strsplit(rawNP,",")[[1]])
        if(!length(np_list)){
          next
        }
        for(j in 1:length(np_list)){
          np_counter<-np_list[[j]]
          #cat(paste("i: ",i," np_counter: ",np_counter," peakArea_num: ",np_df[num_rn,np_counter],"\n"))
          if(is.na(np_df[num_rn,np_list[[j]]])){
            next
          }
          if(is.na(np_df[den_rn,np_list[[j]]])){
            next
          }
          id<-c(id,tempID)
          np<-c(np,1)
          np_ratio=(np_df[num_rn,np_list[[j]]])/(np_df[den_rn,np_list[[j]]])
          ratio<-c(ratio,np_ratio)
        }
      }
      if(length(ratio)!=0){
        answer <- data.frame(id,np,ratio)
        outfile=paste(fn,"/sample",p_list[i],"_timepoint",toString(tp),".txt",sep="")
        write.table(answer,file=outfile,sep="\t",row.names=FALSE,quote=FALSE)
      }
      rm(answer)
      rm(id)
      rm(np)
      rm(ratio)
    }
  }
}

calc_qVals=function(p_df,np_df,link_df,timepts,unpaired,bh){
  
  library(qvalue)
  setwd('C:/Users/knoh1/Documents/qValue_product')
  stats_df<-data.frame()
  fileNames<-character()
  pVals<-numeric()
  qVals<-numeric()
  p_list<-link_df[,2]
  
  for(i in 1:length(p_list)){
    for(j in 1:timepts){
      if(!unpaired){
        #silac stuff here
        tempFileName=paste("silac_df/","sample",colnames(p_df)[i],"_timepoint",j,".txt",sep="")
        if(!file.exists(tempFileName)){
          next
        }
        dat=read.table(tempFileName,sep='\t',header=TRUE)
        model = lm(log(ratio) ~ np, data = dat, na.action = na.omit)
        a=anova(model)
        fileNames<-c(fileNames,tempFileName)
        pVals<-c(pVals,a[['Pr(>F)']][1])
      }else{
        #label free stuff
        tempFileName=paste("labelfree_df/","sample",colnames(p_df)[i],"_timepoint",j,".txt",sep="")
        if(!file.exists(tempFileName)){
          next
        }
        dat=read.table(tempFileName,sep='\t',header=TRUE)
        model1 = lm(log(peakArea) ~ channel + np, data = dat, na.action = na.omit)
        model2 = lm(log(peakArea) ~ channel*np, data = dat, na.action = na.omit)
        a = anova(model1, model2)
        
        fileNames<-c(fileNames,tempFileName)
        pVals<-c(pVals,a[['Pr(>F)']][2])
      }
    }
  }
  
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
  
  lower=0;
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
  
  out_matrix=matrix(qVals,nrow=length(p_list),ncol=4,byrow=TRUE)
  rownames(out_matrix)=p_list
  write.table(out_matrix,file="output.txt",col.names=FALSE)
  
  # phist<-"phist.png"
  # png(filename=phist)
  # hist(pVals,breaks=100)
  # dev.off()
  # qplot<-"qplot.png"
  # png(filename=qplot)
  # qplot=plot(qobj)
  # dev.off()
  # data<-cbind(pVals,qVals)
  
  stats_df <- data.frame(fileNames,pVals,qVals)
  write.table(stats_df,file="master_stats_df.txt",sep="\t",row.names=FALSE,quote=FALSE)
  return(qVals)
}

library(XML)
args = commandArgs(trailingOnly =TRUE)
## The following are arguments that must be passed
p_filePath = "NA"
np_filePath = "NA"
waitPath = "NA"
replicates = 5
timepts = 12
min_rep = 3
## The following are optional arguments that may be passed
label_free = 2
unpaired = 1
silac_num = 1
silac_den = 2
bh = 0
minPA = -Inf
maxPA = Inf
is_min = TRUE
log_LF = FALSE
log_silac = FALSE
is_error = FALSE
p_df<-data.frame()
np_df<-data.frame()

if (length(args) < 5) {
	stop("Invalid Input: Too Few Arguments.", call.=FALSE)
} else if (length(args) > 16) {
	stop("Invalid Input: Too Many Arguments.", call.=FALSE)
} else {
	## The first two arguments when called must be a path to the XML file and a path to the wait file.
	p_filePath = args[1]
	np_filePath = args[2]
	waitPath = args[3]
	## This for loop goes through the rest of the arguments and splits each element into a 1x2 matrix where the elements are
	## the sides of the equals sign (ex. arg_mtx = [timepoints, 12]).
	## Note that for arguments without an equals sign each element in the matrix will be the same (ex. arg_mtx = [min, min]).
	for (i in 4:length(args)){
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
		} else {
			stop("Invalid Input: Unknown Argument.", call.=FALSE)
		}
	}
	p_df<-parseFileHelper(p_filePath, waitPath, replicates, timepts, label_free, unpaired, silac_num, silac_den, min_rep, bh, minPA, maxPA, is_min, log_LF, log_silac, is_error)
	np_df<-parseFileHelper(np_filePath, waitPath, replicates, timepts, label_free, unpaired, silac_num, silac_den, min_rep, bh, minPA, maxPA, is_min, log_LF, log_silac, is_error)
}

#setwd('C:/Users/knoh1/Documents/qValue_ken/offline')
#p_df = read.table('phospho_df.txt',sep=',',header=TRUE,check.names=FALSE)
#np_df = read.table('unphospho_df.txt',sep=',',header=TRUE,check.names=FALSE)
link_df<-parseTextFile()
if(unpaired==0){
  create_silac_df(p_df,np_df,link_df,replicates,timepts,min_rep,label_free,unpaired,silac_num,silac_den,bh,minPA,maxPA,is_min,log_LF,log_silac,is_error)
} else{
  create_lf_df(p_df,np_df,link_df,replicates,timepts,min_rep,label_free,unpaired,silac_num,silac_den,bh,minPA,maxPA,is_min,log_LF,log_silac,is_error)
}
calc_qVals(p_df,np_df,link_df,timepts,unpaired,bh)
