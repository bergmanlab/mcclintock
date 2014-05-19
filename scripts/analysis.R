require(tools)
library(xtable)

bed.files<-list.files(Sys.glob("../sacCer2/*/results"),full.names=T)

for(i in 1:length(bed.files))
{
	if (!exists("dataset"))
	{
		dataset <- read.table(bed.files[i],header=FALSE,col.names=c("Chromosome","Start","End", "TE Name", "Score", "Orientation"),skip=1,sep="\t")
		
		filename<-basename(file_path_sans_ext(bed.files[i]))
		dataset$Sample<-rep(substr(filename,1,9),nrow(dataset))
		dataset$Method<-rep(substr(filename,11,nchar(filename)),nrow(dataset))
	}
	if (exists("dataset"))
	{		
		temp.dataset <- read.table(bed.files[i],header=FALSE,col.names=c("Chromosome","Start","End", "TE Name", "Score", "Orientation"),skip=1,sep="\t")
		
		filename<-basename(file_path_sans_ext(bed.files[i]))
		temp.dataset$Sample<-rep(substr(filename,1,9),nrow(temp.dataset))
		temp.dataset$Method<-rep(substr(filename,11,nchar(filename)),nrow(temp.dataset))
		
		dataset<-rbind(dataset, temp.dataset)
		rm(temp.dataset)
	}
}

new.TEs<-which(dataset$TE.Name == "TY1_new" | dataset$TE.Name == "TY2_new" | dataset$TE.Name == "TY3_new" | dataset$TE.Name == "TY3_1p_new" | dataset$TE.Name == "TY4_new" | dataset$TE.Name == "TY5_new" )

newTEs<-paste(dataset$Method[new.TEs],dataset$Sample[new.TEs])

new.TE.counts<-table(newTEs)

new.TE.counts.ngs_te_mapper<-new.TE.counts[grep("ngs_te_mapper",names(new.TE.counts))]
new.TE.counts.relocate<-new.TE.counts[grep("relocate",names(new.TE.counts))]
new.TE.counts.retroseq<-new.TE.counts[grep("retroseq",names(new.TE.counts))]
new.TE.counts.popoolationte<-new.TE.counts[grep("popoolationte",names(new.TE.counts))]
new.TE.counts.telocate<-new.TE.counts[grep("telocate",names(new.TE.counts))]

retroseq.new.density<-density(new.TE.counts.retroseq)
telocate.new.density<-density(new.TE.counts.telocate)
relocate.new.density<-density(new.TE.counts.relocate)
ngs_te_mapper.new.density<-density(new.TE.counts.ngs_te_mapper)
popoolationte.new.density<-density(new.TE.counts.popoolationte)

pdf('new.counts.pdf')

plot(ngs_te_mapper.new.density,xlim=c(0,750))
lines(telocate.new.density, col = "blue")
lines(relocate.new.density, col = "red")
lines(retroseq.new.density, col = "green")
lines(popoolationte.new.density, col = "purple")

dev.off()

old.TEs<-which(dataset$TE.Name == "TY1_old" | dataset$TE.Name == "TY2_old" | dataset$TE.Name == "TY3_old" | dataset$TE.Name == "TY3_1p_old" | dataset$TE.Name == "TY4_old" | dataset$TE.Name == "TY5_old" )

oldTEs<-paste(dataset$Method[old.TEs],dataset$Sample[old.TEs])

old.TE.counts<-table(oldTEs)

old.TE.counts.ngs_te_mapper<-old.TE.counts[grep("ngs_te_mapper",names(old.TE.counts))]
old.TE.counts.relocate<-old.TE.counts[grep("relocate",names(old.TE.counts))]
old.TE.counts.popoolationte<-old.TE.counts[grep("popoolationte",names(old.TE.counts))]
old.TE.counts.telocate<-old.TE.counts[grep("telocate",names(old.TE.counts))]

telocate.old.density<-density(old.TE.counts.telocate)
relocate.old.density<-density(old.TE.counts.relocate)
ngs_te_mapper.old.density<-density(old.TE.counts.ngs_te_mapper)
popoolationte.old.density<-density(old.TE.counts.popoolationte)

pdf('old.counts.pdf')

plot(ngs_te_mapper.old.density,xlim=c(0,750))
lines(telocate.old.density, col = "blue")
lines(relocate.old.density, col = "red")
lines(popoolationte.old.density, col = "purple")

dev.off()

system("perl concordancecalculation.pl ../sacCer2/*")

new.concordance.files<-list.files(Sys.glob("../sacCer2/*/results"),full.names=T, pattern="new.concordance.tsv$")

for(i in 1:length(new.concordance.files))
{
	if (!exists("new.concordance.dataset"))
	{		
		new.concordance.dataset <- read.table(new.concordance.files[i],header=FALSE,col.names=c("Sample","Method.1","Method.2","Concordance", "High", "Low"),sep="\t")		
	}
	if (exists("new.concordance.dataset"))
	{		
		temp.concordance.dataset <- read.table(new.concordance.files[i],header=FALSE,col.names=c("Sample","Method.1","Method.2","Concordance", "High", "Low"),sep="\t")
		new.concordance.dataset<-rbind(new.concordance.dataset, temp.concordance.dataset)
		rm(temp.concordance.dataset)
	}
}

old.concordance.files<-list.files(Sys.glob("../sacCer2/*/results"),full.names=T, pattern="old.concordance.tsv$")

for(i in 1:length(old.concordance.files))
{
	if (!exists("old.concordance.dataset"))
	{		
		old.concordance.dataset <- read.table(old.concordance.files[i],header=FALSE,col.names=c("Sample","Method.1","Method.2","Concordance", "High", "Low"),sep="\t")	
	}
	if (exists("old.concordance.dataset"))
	{
		temp.concordance.dataset <- read.table(old.concordance.files[i],header=FALSE,col.names=c("Sample","Method.1","Method.2","Concordance", "High", "Low"),sep="\t")
		old.concordance.dataset<-rbind(old.concordance.dataset, temp.concordance.dataset)
		rm(temp.concordance.dataset)
	}
}
	
new.concordance.dataset$Concordance.high <- new.concordance.dataset$Concordance/new.concordance.dataset$High
new.concordance.dataset$Concordance.low <- new.concordance.dataset$Concordance/new.concordance.dataset$Low

new.concordance.mean<-aggregate(cbind(Concordance.high, Concordance.low) ~ Method.1 + Method.2, FUN = mean, data=new.concordance.dataset)
new.concordance.sd<-aggregate(cbind(Concordance.high, Concordance.low) ~ Method.1 + Method.2, FUN = sd, data=new.concordance.dataset)

new.concordance<-cbind(new.concordance.mean, new.concordance.sd$Concordance.high, new.concordance.sd$Concordance.low)
colnames(new.concordance)<-c("Method.1", "Method.2", "Mean.high", "Mean.low", "SD.high", "SD.low")

new.concordance$Method.1<-gsub(" ", "", new.concordance$Method.1)
new.concordance$Method.2<-gsub(" ", "", new.concordance$Method.2)

new.concordance.matrix<-matrix(data = NA, ncol = 5, nrow = 5, dimnames= list(names(table(new.concordance$Method.1)), names(table(new.concordance$Method.1))))
new.methods<-names(table(new.concordance$Method.1))

for (i in 1:length(new.methods))
{
    select<-which(new.concordance$Method.1 == new.methods[i])
    select2<-match(new.concordance$Method.2[select], new.methods)
    new.concordance.matrix[i:length(new.methods),i]<-new.concordance$SD.high[select][select2-i+1]
}

for (i in 1:length(new.methods))
{
    select<-which(new.concordance$Method.1 == new.methods[i])
    select2<-match(new.concordance$Method.2[select], new.methods)
    new.concordance.matrix[i, i:length(new.methods)]<-new.concordance$Mean.high[select][select2-i+1]
}

old.concordance.dataset$Concordance.high <- old.concordance.dataset$Concordance/old.concordance.dataset$High
old.concordance.dataset$Concordance.low <- old.concordance.dataset$Concordance/old.concordance.dataset$Low

old.concordance.mean<-aggregate(cbind(Concordance.high, Concordance.low) ~ Method.1 + Method.2, FUN = mean, data=old.concordance.dataset)
old.concordance.sd<-aggregate(cbind(Concordance.high, Concordance.low) ~ Method.1 + Method.2, FUN = sd, data=old.concordance.dataset)

old.concordance<-cbind(old.concordance.mean, old.concordance.sd$Concordance.high, old.concordance.sd$Concordance.low)
colnames(old.concordance)<-c("Method.1", "Method.2", "Mean.high", "Mean.low", "SD.high", "SD.low")

old.concordance$Method.1<-gsub(" ", "", old.concordance$Method.1)
old.concordance$Method.2<-gsub(" ", "", old.concordance$Method.2)

old.concordance.matrix<-matrix(data = NA, ncol = 4, nrow = 4, dimnames= list(names(table(old.concordance$Method.1)), names(table(old.concordance$Method.1))))
old.methods<-names(table(old.concordance$Method.1))

for (i in 1:length(old.methods))
{
    select<-which(old.concordance$Method.1 == old.methods[i])
    select2<-match(old.concordance$Method.2[select], old.methods)
    old.concordance.matrix[i:length(old.methods),i]<-old.concordance$SD.high[select][select2-i+1]
}

for (i in 1:length(old.methods))
{
    select<-which(old.concordance$Method.1 == old.methods[i])
    select2<-match(old.concordance$Method.2[select], old.methods)
    old.concordance.matrix[i, i:length(old.methods)]<-old.concordance$Mean.high[select][select2-i+1]
}

old.concordance.matrix.table<-xtable(old.concordance.matrix)
new.concordance.matrix.table<-xtable(new.concordance.matrix)

print(old.concordance.matrix.table, floating=FALSE, file="old.concordance.tex")
print(new.concordance.matrix.table, floating=FALSE, file="new.concordance.tex")


