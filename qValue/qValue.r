parseTextFile = function(){
  setwd('C:/Users/knoh1/Documents/qValue_ken')
  link_table = read.table('p_to_np.txt',sep='\t',header=TRUE)
  link_df=data.frame(link_table)
  write.table(link_df,file="link_df.txt",sep="\t")
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
	
	## Print colnames and rownames to check for errors.
	##print(paste0("pXML_df colnames and rownames."))
	##print(paste0(colnames(pXML_df)))
	##print(paste0(rownames(pXML_df)))
	## Creates a copy of the parsed XML data frame to check for errors.
	## write.table(pXML_df, "result_df.csv", row.names=FALSE, col.names=FALSE, sep=",")
	
	## TO DO : LF p Values for both channels.
	ch2pLF = pValueLF(pXML_df, 2, timepts, replicates, min_rep)
	ch3pLF = pValueLF(pXML_df, 3, timepts, replicates, min_rep)
	## TO DO : SILAC p Values. (Note: This call will need to preform a check is SILAC is requested first).
	if (silac_num != -1 && silac_den != -1){
		pSILAC = pValueSILAC(pXML_df, timepts, replicates, min_rep, silac_num, silac_den)
	}
	## TO DO : q Values for LF and SILAC. (Note: The SILAC qValues need to check if SILAC is requested).
	ch2qLF = qValue(ch2pLF)
	ch3qLF = qValue(ch3pLF)
	if (silac_num != -1 && silac_den != -1){
		qSILAC = qValue(pSILAC)
	}
	
	if(identical("C:\\Users\\knoh1\\Documents\\qValue_ken\\phospho_data.xml",filePath)){
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

## TO DO
pValueLF = function(data_df, ch, timepts, replicates, min_rep) {
	
	tpts = timepts/2
	
	## Initialize empty data frame.
	pLF_df = data.frame(matrix(nrow=(tpts),ncol=(dim(data_df)[2])))
	colnames(pLF_df) = colnames(data_df)
	for (i in 1:(tpts)) {
		rownames(pLF_df)[i] = paste("timepoint", toString(i), sep = "")
	}
	
	## Populate data frame with p-values
	apply(data_df, 2, function(col_, ch, tpts, replicates, min_rep){
		
		vec_mtx = matrix(nrow=0, ncol=1)
		low_tp = -1
		low_avg = -Inf
		
		## Create Vectors of timpepoints for p-values.
		for (tp in 1:tpts) {
			vec = NULL
			for (rp in 1:replicates){
				row_name = paste("^", "ch", toString(ch), " rep", toString(rp), " timepoint", toString(tp), "$", sep = "")
				rep_ = as.numeric(col_[grep(paste(row_name, sep = ""), rownames(data_df))])
				## Add data to vector (exclude NA and NaN).
				if(!is.na(rep_)){
					vec = c(vec, rep_)
				}
			}
			
			vec_mtx = t(cbind(t(vec_mtx), matrix(list(vec))))
			
			## Check for lowest timepoint (exclude timepoints without minimum number of replicates.
			if(length(vec) >= min_rep){
				
				
			}
		}
		
		## Do T- Test (Assert that at lest one timepoint has the minimum number of replicates).
		if (low_tp != -1){
			
		}
	}, ch, tpts, replicates, min_rep)
	
	return(pLF_df)
}

## TO DO
pValueSILAC = function(data_df, timepts, replicates, min_rep, silac_num, silac_den) {
	
	tpts = timepts/2
	
	## Initialize empty data frame.
	pSILAC_df = data.frame(matrix(nrow=(tpts),ncol=(dim(data_df)[2])))
	colnames(pSILAC_df) = colnames(data_df)
	for (i in 1:(tpts)) {
		rownames(pSILAC_df)[i] = paste("timepoint", toString(i), sep = "")
	}
	
	## Print colnames and rownames to check for errors.
	##print(paste0("pValueSILAC colnames and rownames."))
	##print(paste0(colnames(pSILAC_df)))
	##print(paste0(rownames(pSILAC_df)))
	
	return(pSILAC_df)
}

## TO DO
qValue = function(pVal_df) {
	
	## Initialize empty data frame.
	qVal_df = data.frame(matrix(nrow=(dim(pVal_df)[1]),ncol=(dim(pVal_df)[2])))
	colnames(qVal_df) = colnames(pVal_df)
	rownames(qVal_df) = rownames(pVal_df)
	
	## Print colnames and rownames to check for errors.
	##print(paste0("qValue colnames and rownames."))
	##print(paste0(colnames(qVal_df)))
	##print(paste0(rownames(qVal_df)))
	
	return(qVal_df)
}

library(XML)
args = commandArgs(trailingOnly =TRUE)
## The following are arguments that must be passed
p_filePath = "NA"
np_filePath = "NA"
waitPath = "NA"
replicates = 5
timepts = 16
min_rep = 3
## The following are optional arguments that may be passed
label_free = 1
unpaired = -1
silac_num = -1
silac_den = -1
bh = 0
minPA = -Inf
maxPA = Inf
is_min = TRUE
log_LF = FALSE
log_silac = FALSE
is_error = FALSE

cat(args,file="args.txt",sep="\n")

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
	link_df<-parseTextFile()
	p_df<-parseFileHelper(p_filePath, waitPath, replicates, timepts, label_free, unpaired, silac_num, silac_den, min_rep, bh, minPA, maxPA, is_min, log_LF, log_silac, is_error)
	np_df<-parseFileHelper(np_filePath, waitPath, replicates, timepts, label_free, unpaired, silac_num, silac_den, min_rep, bh, minPA, maxPA, is_min, log_LF, log_silac, is_error)
}
