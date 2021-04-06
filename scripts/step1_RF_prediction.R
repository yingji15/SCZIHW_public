
#use: Rscript ./rf_0212.R /fs0/jiy1/autism/features/database/expr_network/Rscript/depict/depict_hgncsym.rds  /fs0/chenr6/gibbs/file/scz.top1.v7.t0 /fs0/chenr6/gibbs/file/LPG.v7.t0 2




args<-commandArgs(TRUE)
if (length(args)==0){
  stop("at least one argument must be supplied.n",call. = FALSE)
}


#####load data#####
fea.file<-args[1]



#here it should be: "/fs0/jiy1/autism/features/database/expr_network/Rscript/depict"
mylist<-unlist(strsplit(fea.file,"[/]") )
feaname<-mylist[length(mylist)]
print(feaname)
pos.file<-args[2]
neg.file<-args[3]

print(pos.file)
#save positive name for naming
posname<-sapply( strsplit(pos.file,"/"),'[',length(unlist(strsplit(pos.file,"/"))) )

print(strsplit(posname,"\\."))

pos2<-sapply( strsplit(posname,"\\."), '[',1)
print(pos2)
#note down write to path for prediction file
#concept: use sub by matching a / followed by zero or more characters that are not a / till the end ($) of the string and replace it with blank ('')
#writepath<-paste(  paste( sub('/[^/]*$','',fea.file), pos2,sep="_"),"/pred/",sep="")

#if want to save to feature file directory
#writepath<-paste(  paste( sub('/[^/]*$','',fea.file),"/pred/",sep="") )

#if want to save to pos file directory
#writepath<-paste(  paste( sub('/[^/]*$','',pos.file),"/pred/",sep="") )
cwd<-getwd()
writepath<-paste(cwd,"/pred/",sep="")
print(writepath)

#check if the write path exist and create if not

if (!dir.exists(writepath)){
dir.create(writepath)
} else {
    print("output dir already exists!")
}

print(writepath)
setwd(writepath)
#how many times you want to run it
ntimes<-as.numeric(args[4])

print(fea.file)
print(pos.file)
print(neg.file)
print(ntimes)

filetype<-unlist(strsplit(feaname,"[.]"))[-1]
print(filetype)
	if (filetype==c("rds")) {
#for depict data, used rds
	pre.fea<-readRDS(fea.file)
	} else {
#for other data, sue read.table
	pre.fea<-read.table(fea.file,head=T,as.is=T,fill=T,row.names=NULL)
	}
	
colnames(pre.fea)[1]<-c("gene")

fea<-pre.fea

pos<-read.table(pos.file,head=F,as.is=T,fill=T)
#pos<-read.table("/fs0/chenr6/gibbs/file/scz.top1.v7.t0",head=F,as.is=T,fill=T)
neg<-read.table(neg.file,head=F,as.is=T,fill=T)
#neg<-read.table("/fs0/chenr6/gibbs/file/LPG.v7.t0",as.is=T,head=F,fill=T)


#print(dim(fea))
#print(dim(pos))


posfea<-fea[fea[,1] %in% pos[,1],]
negfea<-fea[fea[,1] %in% neg[,1],]
#change here: gen fea is all now
genfea<-fea
genfea[is.na(genfea)]<-0
genname<-as.character(genfea$gene)
print("gene name")
print(genname[1:5])
print(dim(genfea))

#combine positive and negative set
setfea<-rbind(posfea,negfea)
#reset rownames (so index will be in order)
#rownames(setfea)<-NULL
#setnames<-setfea[,1]
#print(setnames)
#print(colnames(setfea)[1:5])

#output
genpred<-list()
rfresult<-list()


setlabels<-c(rep(1,nrow(posfea)),rep(0,nrow(negfea)))

library(randomForest)
library("pROC")


rftestroc<-svmtestroc<-enettestroc<-rep(0,ntimes)
out<-data.frame() #store misclassification results
genout<-data.frame() #store genome wide predictions

for (i in 1:ntimes){
##set up
	
	##downsample negative (make pos and neg equal)
	selectpos<-sample(1:nrow(posfea),nrow(posfea))
	selectneg<-sample(1:nrow(negfea),nrow(posfea))
	
	x<-rbind(posfea[selectpos,],negfea[selectneg,])
	rownames(x)<-NULL #reset row index
	xgene<-x[,1] #store names
	x<-x[,-1] #delete name column
	
	#construct y label
	y<-factor(c(rep(1,nrow(posfea)),rep(0,nrow(posfea))))
	##test,train split
	train <- sample(1:nrow(x), nrow(x) * 2/3);
	test <- (-train);
	print(dim(x))
	
	
	md.rf <- randomForest(x[train,], y[train], importance=TRUE, proximity=TRUE,ntree=3000)
	results.rf<-predict(md.rf,x,type="prob")
	
	rf.train_roc<-roc(y[train],results.rf[,2][train],direction="<");
	rf.test_roc <- roc(y[test], results.rf[,2][test],direction="<");
	rf.all_roc <- roc(y, results.rf[,2],direction="<");
	
	if (i==1){
		plotfile=paste(writepath,feaname,"_auc.pdf",sep="")
		pdf(file=plotfile)
		par(las = 2, mai = c(1.5, 1.3, 1.0, 0.10))
	# plot(train_roc, col="blue", cex.axis=2, lwd.ticks=2,cex.lab=1.5, lwd =5)
		plot(rf.test_roc,  print.auc=TRUE, col="red", cex.axis=2, lwd.ticks=2,cex.lab=1.5, lwd =5);
	# plot(all_roc, add=TRUE, cex.axis=2, lwd.ticks=2,cex.lab=1.5, lwd =5);
		legend("right",cex=2,c("test"),title="ROC Curve",lwd=5,col=c("red"),inset = .02);
		dev.off()
		}
	

	#plot several feature imp	
	#rfoptions <- round(importance(md.rf), 2)[which(round(importance(md.rf), 2)[,2] != 0),]
	if (i <=3){
	pdf(file=paste(writepath,feaname,i,"_allfeaimp.pdf",sep=""))
	varImpPlot(md.rf,type=2)
	dev.off()
#	write.table(rfoptions,file = paste("rf.options",i,sep="."));
	#rfresult[[i]]<-round( importance(md.rf), 2)	
	}


##save feature importance
	feanames<-rownames(importance(md.rf))
	rfresult[[i]]<-round( importance(md.rf,type=2),2)
	rftestroc[i]<-as.numeric(pROC::auc(rf.test_roc))
	
	
	
 

##genome wide prediction
	gen.rf <- randomForest(x, y, importance=TRUE, proximity=TRUE,ntree=3000)
	genresults.rf<-predict(gen.rf,genfea,type="prob")[,2]
		
	genpred[[i]]<-genresults.rf
	print(length(genresults.rf))
	
	}	
	
	#genome pred
	bigpred=do.call(cbind,genpred)
	final=cbind(genname,bigpred)

	#combine feature importance
	bigrf=do.call(cbind,rfresult)
	rffinal=cbind(feanames,bigrf) 
	
	
print("auc roc")
print(rftestroc)
print(mean(rftestroc))

###write to a file 
write.table(final,file=paste(writepath,feaname,pos2,ntimes,"allgen_pred.txt", sep=""),row.names=FALSE,quote=F)
#write.table(out,file="/fs0/jiy1/autism/features/database/go_network/gtex5050/summ_wrong_pos.txt",quote=F)

#write.table(rftestroc,file=paste(writepath,feaname,pos2,ntimes, "alltestauc.txt",sep=""),row.names=FALSE,quote=F)

##write importance to a file
#write.table(rffinal,file=paste(writepath,paste( feaname,pos2,ntimes,"allcomb_fea_imp.txt",sep="_"),sep=""),quote=F,row.names=FALSE)
