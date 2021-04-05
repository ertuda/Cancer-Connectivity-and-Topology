# Analysis of the Co-expression Networks of 
# Colorectal Adenoma and Carcinoma Datasets
# By Ertugrul Dalgic, PhD
# 2015-2018

# number of randomizations
norands = 1000
# randomization p value threshold
randpthr = 0.05

# index of normal and tumor samples
e8671ns = 1:32
e8671ts = 33:64
e18105ns = 1:17
e18105ts = 18:34

# Normalize the raw datasets with RMA from affy library
library(affy)

# Colorectal Adenoma Dataset (GSE8671)
setwd('C:/Users/Dell/Documents/Projects/Colorectal Topology/GSE8671/')
eade <- read.AnnotatedDataFrame('GSE8671_targets.txt',header=FALSE,
                                row.names=1,as.is=TRUE)
raweade = ReadAffy(filenames=pData(eade)$FileName,phenoData=eade)
normeade = rma(raweade)
e8671 = exprs(normeade)
remove(normeade);remove(raweade)

# Colorectal Carcinoma Dataset (GSE18105)
setwd('C:/Users/Dell/Documents/Projects/Colorectal Topology/GSE18105/')
ecar <- read.AnnotatedDataFrame('GSE18105_targets.txt',header=FALSE,
                                row.names=1,as.is=TRUE)
rawecar = ReadAffy(filenames=pData(ecar)$FileName,phenoData=ecar)
normecar = rma(rawecar)
e18105 = exprs(normecar)
remove(normecar);remove(rawecar)

setwd('C:/Users/Dell/Documents/Projects/Colorectal Topology/')

# first reducing affy ids
# collapse by max var
# GSE8671 and GSE18105 are both Affymetrix Hgu133 Plus2
library(hgu133plus2.db)
mapped_probes <- mappedkeys(hgu133plus2ENTREZID)
hid <- as.list(hgu133plus2ENTREZID[mapped_probes])

e8671gid <- e8671[intersect(rownames(e8671),names(hid)),]
rownames(e8671gid) <- sapply(rownames(e8671gid),function(i) hid[[i]])
# use var to select the probe set
# use max of normal or tumor var
axn <- apply(e8671gid[,e8671ns],1,sd)
axt <- apply(e8671gid[,e8671ts],1,sd)
ax <- apply(cbind(axn,axt),1,max)
sax <- split(ax,as.factor(names(ax)))
ssax <- sapply(sax,max)
sind <- sapply(seq_along(ax),function(i) ax[i]==ssax[names(ax)[i]])
e8671j <- e8671gid[sind,]
remove(e8671gid)
remove(e8671)

e18105gid <- e18105[intersect(rownames(e18105),names(hid)),]
rownames(e18105gid) <- sapply(rownames(e18105gid),function(i) hid[[i]])
# use var to select the probe set
axn <- apply(e18105gid[,e18105ns],1,sd)
axt <- apply(e18105gid[,e18105ts],1,sd)
ax <- apply(cbind(axn,axt),1,max)
sax <- split(ax,as.factor(names(ax)))
ssax <- sapply(sax,max)
sind <- sapply(seq_along(ax),function(i) ax[i]==ssax[names(ax)[i]])
e18105j <- e18105gid[sind,]
remove(e18105gid)
remove(e18105)

# next, we reduce genes by average expression level
# to avoid lowly expressed genes from correlation based analysis
lowexpthr = 5
me8671j <- rowMeans(e8671j)
me18105j <- rowMeans(e18105j)
e8671jh <- e8671j[me8671j>=lowexpthr,]
e18105jh <- e18105j[me18105j>=lowexpthr,]

# next, we reduce genes by variance
# to avoid lowly varied genes from correlation based analysis
lowvarthr = 0.5
de8671jhn <- apply(e8671jh[,e8671ns],1,sd)
de8671jht <- apply(e8671jh[,e8671ts],1,sd)
de8671jh <- apply(cbind(de8671jhn,de8671jht),1,max)
de18105jhn <- apply(e18105jh[,e18105ns],1,sd)
de18105jht <- apply(e18105jh[,e18105ts],1,sd)
de18105jh <- apply(cbind(de18105jhn,de18105jht),1,max)
e8671jhv <- e8671jh[de8671jh>=lowvarthr,]
e18105jhv <- e18105jh[de18105jh>=lowvarthr,]
print(paste(nrow(e8671jhv),' genes are selected for GSE8671'))
print(paste(nrow(e18105jhv),' genes are selected for GSE18105'))

remove(e8671j);remove(e18105j)
remove(e8671jh);remove(e18105jh)

# finally, calculate spearman correlations and 
# multiple testing corrected p-values
library(Hmisc)
library(metaSEM)

# Calculation of weighted degree values for the cancer networks
# self correlations are set to 0
# positive and negative correlations are considered separately
# connection requires significant correlation based on randomization
lowthr = 0.2
highthr = 0.5

# real correlations based on real normal and tumor sample sets
# for e8671 (adenoma), samples (columns) 1-32  are normal, 33-64 are tumor
# for e18105 (carcinoma 1), samples (columns) 1-17 are normal, 18-34 are tumor
noradeexp <- e8671jhv[,e8671ns];tumadeexp <- e8671jhv[,e8671ts]
norcarexp <- e18105jhv[,e18105ns]; tumcarexp <- e18105jhv[,e18105ts]

# adenoma
ncor <- rcorr(t(noradeexp),type='pearson')

aposn1 <- ncor$r
diag(aposn1) = 0
aposn1[aposn1 < lowthr] = 0
aposn1[aposn1 != 0] = 1
anegn1 <- ncor$r
diag(anegn1) = 0
anegn1[anegn1 > -lowthr] = 0
anegn1[anegn1 != 0] = 1

aposn2 <- ncor$r
diag(aposn2) = 0
aposn2[aposn2 < highthr] = 0
aposn2[aposn2 != 0] = 1
anegn2 <- ncor$r
diag(anegn2) = 0
anegn2[anegn2 > -highthr] = 0
anegn2[anegn2 != 0] = 1

# to be used for randoms later
raposn <- ncor$r
raposn[] = 0
ranegn <- ncor$r
ranegn[] = 0

tcor <- rcorr(t(tumadeexp),type='pearson')

apost1 <- tcor$r
diag(apost1) = 0
apost1[apost1 < lowthr] = 0
apost1[apost1 != 0] = 1
anegt1 <- tcor$r
diag(anegt1) = 0
anegt1[anegt1 > -lowthr] = 0
anegt1[anegt1 != 0] = 1

apost2 <- tcor$r
diag(apost2) = 0
apost2[apost2 < highthr] = 0
apost2[apost2 != 0] = 1
anegt2 <- tcor$r
diag(anegt2) = 0
anegt2[anegt2 > -highthr] = 0
anegt2[anegt2 != 0] = 1

# to be used for randoms later
rapost <- tcor$r
rapost[] = 0
ranegt <- tcor$r
ranegt[] = 0

# carcinoma
ncor <- rcorr(t(norcarexp),type='pearson')

cposn1 <- ncor$r
diag(cposn1) = 0
cposn1[cposn1 < lowthr] = 0
cposn1[cposn1 != 0] = 1
cnegn1 <- ncor$r
diag(cnegn1) = 0
cnegn1[cnegn1 > -lowthr] = 0
cnegn1[cnegn1 != 0] = 1

cposn2 <- ncor$r
diag(cposn2) = 0
cposn2[cposn2 < highthr] = 0
cposn2[cposn2 != 0] = 1
cnegn2 <- ncor$r
diag(cnegn2) = 0
cnegn2[cnegn2 > -highthr] = 0
cnegn2[cnegn2 != 0] = 1

# to be used for randoms later
rcposn <- ncor$r
rcposn[] = 0
rcnegn <- ncor$r
rcnegn[] = 0

tcor <- rcorr(t(tumcarexp),type='pearson')

cpost1 <- tcor$r
diag(cpost1) = 0
cpost1[cpost1 < lowthr] = 0
cpost1[cpost1 != 0] = 1
cnegt1 <- tcor$r
diag(cnegt1) = 0
cnegt1[cnegt1 > -lowthr] = 0
cnegt1[cnegt1 != 0] = 1

cpost2 <- tcor$r
diag(cpost2) = 0
cpost2[cpost2 < highthr] = 0
cpost2[cpost2 != 0] = 1
cnegt2 <- tcor$r
diag(cnegt2) = 0
cnegt2[cnegt2 > -highthr] = 0
cnegt2[cnegt2 != 0] = 1

# to be used for randoms later
rcpost <- ncor$r
rcpost[] = 0
rcnegt <- ncor$r
rcnegt[] = 0

# randomized samples based correlations
# take random sample sets of size 32 and 17
for (jj in 1:norands){
  print(jj)
  
  sc = c(1:e8671ts[length(e8671ts)])
  randset = sample(sc,e8671ns[length(e8671ns)])
  randset2 = sc[-randset]
  randnoradeexp <- e8671jhv[,randset];randtumadeexp <- e8671jhv[,randset2]
  
  sc = c(1:e18105ts[length(e18105ts)])
  randset = sample(sc,e18105ns[length(e18105ns)])
  randset2 = sc[-randset]
  randnorcarexp <- e18105jhv[,randset];randtumcarexp <- e18105jhv[,randset2]
  
  # adenoma
  ncor <- rcorr(t(randnoradeexp),type='pearson')
  
  aposn <- ncor$r
  diag(aposn) = 0
  aposn[aposn < highthr] = 0
  aposn[aposn != 0] = 1
  anegn <- ncor$r
  diag(anegn) = 0
  anegn[anegn > -highthr] = 0
  anegn[anegn != 0] = 1
  raposn = raposn + aposn
  ranegn = ranegn + anegn
  
  tcor <- rcorr(t(randtumadeexp),type='pearson')

  apost <- tcor$r
  diag(apost) = 0
  apost[apost < highthr] = 0
  apost[apost != 0] = 1
  anegt <- tcor$r
  diag(anegt) = 0
  anegt[anegt > -highthr] = 0
  anegt[anegt != 0] = 1
  rapost = rapost + apost
  ranegt = ranegt + anegt
  
  # carcinoma
  ncor <- rcorr(t(randnorcarexp),type='pearson')

  cposn <- ncor$r
  diag(cposn) = 0
  cposn[cposn < highthr] = 0
  cposn[cposn != 0] = 1
  cnegn <- ncor$r
  diag(cnegn) = 0
  cnegn[cnegn > -highthr] = 0
  cnegn[cnegn != 0] = 1
  rcposn = rcposn + cposn
  rcnegn = rcnegn + cnegn
  
  tcor <- rcorr(t(randtumcarexp),type='pearson')
  
  cpost <- tcor$r
  diag(cpost) = 0
  cpost[cpost < highthr] = 0
  cpost[cpost != 0] = 1
  cnegt <- tcor$r
  diag(cnegt) = 0
  cnegt[cnegt > -highthr] = 0
  cnegt[cnegt != 0] = 1
  rcpost = rcpost + cpost
  rcnegt = rcnegt + cnegt
}

# pvalues based on randomizations
# FDR correction of randomization based p values
raposn = raposn/norands
rax = as.vector(raposn)[as.vector(upper.tri(raposn))]
x = p.adjust(rax,method='fdr')
raposn = vec2symMat(x, diag = FALSE, byrow = TRUE)
diag(raposn) = 0
raposn[raposn < randpthr] = 0
raposn[raposn != 0] = 1

rapost = rapost/norands
rax = as.vector(rapost)[as.vector(upper.tri(rapost))]
x = p.adjust(rax,method='fdr')
rapost = vec2symMat(x, diag = FALSE, byrow = TRUE)
diag(rapost) = 0
rapost[rapost < randpthr] = 0
rapost[rapost != 0] = 1

ranegn = ranegn/norands
rax = as.vector(ranegn)[as.vector(upper.tri(ranegn))]
x = p.adjust(rax,method='fdr')
ranegn = vec2symMat(x, diag = FALSE, byrow = TRUE)
diag(ranegn) = 0
ranegn[ranegn < randpthr] = 0
ranegn[ranegn != 0] = 1

ranegt = ranegt/norands
rax = as.vector(ranegt)[as.vector(upper.tri(ranegt))]
x = p.adjust(rax,method='fdr')
ranegt = vec2symMat(x, diag = FALSE, byrow = TRUE)
diag(ranegt) = 0
ranegt[ranegt < randpthr] = 0
ranegt[ranegt != 0] = 1

rcposn = rcposn/norands
rax = as.vector(rcposn)[as.vector(upper.tri(rcposn))]
x = p.adjust(rax,method='fdr')
rcposn = vec2symMat(x, diag = FALSE, byrow = TRUE)
diag(rcposn) = 0
rcposn[rcposn < randpthr] = 0
rcposn[rcposn != 0] = 1

rcpost = rcpost/norands
rax = as.vector(rcpost)[as.vector(upper.tri(rcpost))]
x = p.adjust(rax,method='fdr')
rcpost = vec2symMat(x, diag = FALSE, byrow = TRUE)
diag(rcpost) = 0
rcpost[rcpost < randpthr] = 0
rcpost[rcpost != 0] = 1

rcnegn = rcnegn/norands
rax = as.vector(rcnegn)[as.vector(upper.tri(rcnegn))]
x = p.adjust(rax,method='fdr')
rcnegn = vec2symMat(x, diag = FALSE, byrow = TRUE)
diag(rcnegn) = 0
rcnegn[rcnegn < randpthr] = 0
rcnegn[rcnegn != 0] = 1

rcnegt = rcnegt/norands
rax = as.vector(rcnegt)[as.vector(upper.tri(rcnegt))]
x = p.adjust(rax,method='fdr')
rcnegt = vec2symMat(x, diag = FALSE, byrow = TRUE)
diag(rcnegt) = 0
rcnegt[rcnegt < randpthr] = 0
rcnegt[rcnegt != 0] = 1


# raposn, rapost, .... capture random connections so they will be substracted  
# posn2, post2, negn2 and negt2 ignore the correlation on the other sample type
# differential correlation: hardposn, hardpost, hardnegn and hardnegt
# correlation >= highthr in the sample and lower than 0.05 in the randomized correlations
# and <= lowthr in corresponding normal or tumor

aposn2 = aposn2 - raposn
hardposn = aposn2 - apost1
hardposn[hardposn < 0] = 0

apost2 = apost2 - rapost
hardpost = apost2 - aposn1
hardpost[hardpost < 0] = 0

anegn2 = anegn2 - ranegn
hardnegn = anegn2 - anegt1
hardnegn[hardnegn < 0] = 0

anegt2 = anegt2 - ranegt
hardnegt = anegt2 - anegn1
hardnegt[hardnegt < 0] = 0

napd <- rowSums(hardposn);nand <- rowSums(hardnegn)
tapd <- rowSums(hardpost);tand <- rowSums(hardnegt)

cposn2 = cposn2 - rcposn
hardposn = cposn2 - cpost1
hardposn[hardposn < 0] = 0

cpost2 = cpost2 - rcpost
hardpost = cpost2 - cposn1
hardpost[hardpost < 0] = 0

cnegn2 = cnegn2 - rcnegn
hardnegn = cnegn2 - cnegt1
hardnegn[hardnegn < 0] = 0

cnegt2 = cnegt2 - rcnegt
hardnegt = cnegt2 - cnegn1
hardnegt[hardnegt < 0] = 0

ncpd <- rowSums(hardposn);ncnd <- rowSums(hardnegn)
tcpd <- rowSums(hardpost);tcnd <- rowSums(hardnegt)

# less and more connected genes based on their weighted degree difference
difapd <- tapd-napd
lessapd <- difapd[difapd<0];lessapd <- abs(lessapd)
moreapd <- difapd[difapd>0]
difcpd <- tcpd-ncpd
lesscpd <- difcpd[difcpd<0];lesscpd <- abs(lesscpd)
morecpd <- difcpd[difcpd>0]

difand <- tand-nand
lessand <- difand[difand<0];lessand <- abs(lessand)
moreand <- difand[difand>0]
difcnd <- tcnd-ncnd
lesscnd <- difcnd[difcnd<0];lesscnd <- abs(lesscnd)
morecnd <- difcnd[difcnd>0]

# write less and more connected genes with weights based on weighted degree difference
write.table(lessapd,file='lessapd.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(moreapd,file='moreapd.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lesscpd,file='lesscpd.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(morecpd,file='morecpd.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessand,file='lessand.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(moreand,file='moreand.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lesscnd,file='lesscnd.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(morecnd,file='morecnd.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)

# write top ranking genes
lessapd5 = lessapd[lessapd >= quantile(lessapd,probs=0.95)]
lessapd10 = lessapd[lessapd >= quantile(lessapd,probs=0.9)]
lessapd20 = lessapd[lessapd>=quantile(lessapd,probs=0.8)] 
write.table(lessapd5,file='lessapd_top5percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessapd10,file='lessapd_top10percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessapd20,file='lessapd_top20percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)

perc5 = quantile(moreapd,probs=0.95)
perc10 = quantile(moreapd,probs=0.9)
perc20 = quantile(moreapd,probs=0.8)
write.table(moreapd[moreapd>=perc5],file='moreapd_top5percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(moreapd[moreapd>=perc10],file='moreapd_top10percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(moreapd[moreapd>=perc20],file='moreapd_top20percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)

perc5 = quantile(lesscpd,probs=0.95)
perc10 = quantile(lesscpd,probs=0.9)
perc20 = quantile(lesscpd,probs=0.8)
write.table(lesscpd[lesscpd>=perc5],file='lesscpd_top5percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lesscpd[lesscpd>=perc10],file='lesscpd_top10percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lesscpd[lesscpd>=perc20],file='lesscpd_top20percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)

perc5 = quantile(morecpd,probs=0.95)
perc10 = quantile(morecpd,probs=0.9)
perc20 = quantile(morecpd,probs=0.8)
write.table(morecpd[morecpd>=perc5],file='morecpd_top5percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(morecpd[morecpd>=perc10],file='morecpd_top10percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(morecpd[morecpd>=perc20],file='morecpd_top20percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)

perc5 = quantile(lessand,probs=0.95)
perc10 = quantile(lessand,probs=0.9)
perc20 = quantile(lessand,probs=0.8)
write.table(lessand[lessand>=perc5],file='lessand_top5percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessand[lessand>=perc10],file='lessand_top10percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessand[lessand>=perc20],file='lessand_top20percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)

perc5 = quantile(moreand,probs=0.95)
perc10 = quantile(moreand,probs=0.9)
perc20 = quantile(moreand,probs=0.8)
write.table(moreand[moreand>=perc5],file='moreand_top5percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(moreand[moreand>=perc10],file='moreand_top10percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(moreand[moreand>=perc20],file='moreand_top20percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)

perc5 = quantile(lesscnd,probs=0.95)
perc10 = quantile(lesscnd,probs=0.9)
perc20 = quantile(lesscnd,probs=0.8)
write.table(lesscnd[lesscnd>=perc5],file='lesscnd_top5percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lesscnd[lesscnd>=perc10],file='lesscnd_top10percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lesscnd[lesscnd>=perc20],file='lesscnd_top20percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)

perc5 = quantile(morecnd,probs=0.95)
perc10 = quantile(morecnd,probs=0.9)
perc20 = quantile(morecnd,probs=0.8)
write.table(morecnd[morecnd>=perc5],file='morecnd_top5percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(morecnd[morecnd>=perc10],file='morecnd_top10percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(morecnd[morecnd>=perc20],file='morecnd_top20percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)

# neighbors of top ranking genes
# write top ranking genes
perc5 = quantile(lessapd,probs=0.95)
perc10 = quantile(lessapd,probs=0.9)
perc20 = quantile(lessapd,probs=0.8)
write.table(lessapd[lessapd>=perc5],file='lessapd_top5percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessapd[lessapd>=perc10],file='lessapd_top10percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessapd[lessapd>=perc20],file='lessapd_top20percent.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)



# plot degs and deg differences
#png(file='adeno_carci_degs.png',width=480,height=800)
tiff(file='adeno_carci_degs.tiff',width=3.35, height=6,units="in",res=600)
par(mfrow = c(2,2),mar=c(2.3,3.9,2.3,1))

# plot ABSOLUTE degree values
boxplot(napd,tapd,main='Ade 1 Pos',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5), xaxt='n',
        cex.axis=1.15,cex.lab=1.25)
axis(1,at=c(0.6,1.5),labels=c('N','T'),tick=FALSE,cex.axis=1.15)

boxplot(ncpd,tcpd,main='Car 1 Pos',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5), xaxt='n',
        cex.axis=1.15,cex.lab=1.25)
axis(1,at=c(0.6,1.5),labels=c('N','T'),tick=FALSE,cex.axis=1.15)

boxplot(nand,tand,main='Ade 1 Neg',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5), xaxt='n',
        cex.axis=1.15,cex.lab=1.25)
axis(1,at=c(0.6,1.5),labels=c('N','T'),tick=FALSE,cex.axis=1.15)

boxplot(ncnd,tcnd,main='Car 1 Neg',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5), xaxt='n',
        cex.axis=1.15,cex.lab=1.25)
axis(1,at=c(0.6,1.5),labels=c('N','T'),tick=FALSE,cex.axis=1.15)
dev.off()

# mann whitney tests for comparing degree dist. of tumor to normal (paired)
wt <- wilcox.test(napd,tapd,paired=TRUE,alternative='greater')
print(paste('adenoma positive paired wilcox test p-value:',wt$p.value))
wt <- wilcox.test(nand,tand,paired=TRUE,alternative='greater')
print(paste('adenoma negative paired wilcox test p-value:',wt$p.value))
wt <- wilcox.test(ncpd,tcpd,paired=TRUE,alternative='greater')
print(paste('carcinoma positive paired wilcox test p-value:',wt$p.value))
wt <- wilcox.test(ncnd,tcnd,paired=TRUE,alternative='greater')
print(paste('carcinoma negative paired wilcox test p-value:',wt$p.value))

##############################################################
# plot the adenoma and carcinoma expression values of a gene pair of interest 
# adenoma
givengene1 = '6352'
givengene2 = '2049'
par(mfrow=c(1,2))
plot(noradeexp[givengene1,],noradeexp[givengene2,],
     xlab=givengene1,ylab=givengene2,main='adenoma normal',cex.lab=1.5,cex.main=1.5)
plot(tumadeexp[givengene1,],tumadeexp[givengene2,],
     xlab=givengene1,ylab=givengene2,main='adenoma tumor',cex.lab=1.5,cex.main=1.5)
print('normal Pearson correlation test')
cor.test(noradeexp[givengene1,],noradeexp[givengene2,],method='pearson')
print('tumor Pearson correlation test')
cor.test(tumadeexp[givengene1,],tumadeexp[givengene2,],method='pearson')

# carcinoma
givengene1 = '6352'
givengene2 = '2049'
par(mfrow=c(1,2))
plot(norcarexp[givengene1,],norcarexp[givengene2,],
     xlab=givengene1,ylab=givengene2,main='carcinoma normal',cex.lab=1.5,cex.main=1.5)
plot(tumcarexp[givengene1,],tumcarexp[givengene2,],
     xlab=givengene1,ylab=givengene2,main='carcinoma tumor',cex.lab=1.5,cex.main=1.5)
print('normal Pearson correlation test')
cor.test(norcarexp[givengene1,],norcarexp[givengene2,],method='pearson')
print('tumor Pearson correlation test')
cor.test(tumcarexp[givengene1,],tumcarexp[givengene2,],method='pearson')

