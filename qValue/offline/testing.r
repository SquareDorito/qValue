getInput = function(){
  setwd('C:/Users/knoh1/Documents/qValue_ken/offline')
  link_table = read.table('p_to_np.txt',sep='\t',header=TRUE)
  link_df=data.frame(link_table)
  rownames(link_df)<-link_df[,2]
  np<-as.character(link_df[,3][1])
  #rownames(link_df)
  #cat(paste("list of np: ",np,"\n"))
  np_list<-as.list(strsplit(np,",")[[1]])
  #cat(paste(length(np_list),"  ",np,"  ",np_list,"\n"))
  #np_list[1]
  np_list[[3]]
  #cat(paste(np_list[1],"\n",np_list[2],"\n",np_list[3]))
  
  #link_df[[1]][2]
  #colnames(link_df)[1]
  # num_entries<-vector("list",12)
  # num_entries[1]<-1
  # num_entries[3]<-3
  # num_entries[5]<-5
  # 
  # ex<-c(1,2,3,4,5)
  # ex<-c(ex,6)
  # ex

}

testPhospho=function(){
  setwd('C:/Users/knoh1/Documents/qValue_ken/offline')
  link_table = read.table('p_to_np.txt',sep='\t',header=TRUE)
  p_df = read.table('phospho_df.txt',sep=',',header=TRUE,check.names=FALSE)
  link_df=data.frame(link_table)
  rownames(link_df)<-link_df[,2]
  np<-as.character(link_df[,3][1])
  np_list<-as.list(strsplit(np,",")[[1]])
  np_list[[3]]<-32338205
  cat(" ",paste(p_df['ch1 rep1 timepoint1',np_list[[3]]]))
 
  # if(is.na(p_df['ch1 rep1 timepoint1',np_list[[3]]])){
  #   cat("hello")
  # }
  #p_df["ch1 rep1 timepoint1",c("X32338205")]
}

testPhospho()