# Analysis of the Co-expression Networks of 
# Aldosterone producing adenoma  Dataset
# By Ertugrul Dalgic, PhD
# 2019

# number of randomizations
norands = 1000
# randomization p value threshold
randpthr = 0.05

# index of normal and tumor samples
e60042adns = 1:7
e60042adts = 8:14

# Already log transformed and quantile normalized
# first colon is agilent id: will be converted to entrez gene id
# Before here, ID columnd is Agi IDs
# After this section ID column is GeneIDs
library(preprocessCore)
setwd('C:/Users/Dell/Documents/Projects/Colorectal Topology/GSE60042/')
agd <- read.table('GPL14550-GeneID.txt',header=TRUE)

e60042ad <- read.table('GSE60042_adenoma.txt',header=TRUE,as.is=TRUE)
e60042ad$ID <- sapply(e60042ad$ID,function(i) agd[agd$AgiID==i,]$GeneID)

# first reducing ids
# collapse by max var
# use var to select the probe set
# use max of normal or tumor var
# First colon is IDs, which is moved as rowname at last
setwd('C:/Users/Dell/Documents/Projects/Colorectal Topology/')

axn <- apply(e60042ad[,e60042adns+1],1,sd)
axt <- apply(e60042ad[,e60042adts+1],1,sd)
ax <- apply(cbind(axn,axt),1,max)
names(ax) = e60042ad$ID
sax <- split(ax,as.factor(names(ax)))
ssax <- sapply(sax,max)
sind <- sapply(seq_along(ax),function(i) ax[i]==ssax[names(ax)[i]])
e60042ad <- e60042ad[sind,]

rownames(e60042ad) = e60042ad$ID
e60042ad <- e60042ad[,-1]

# next, we reduce genes by average expression level
# to avoid lowly expressed genes from correlation based analysis
lowexpthr = 5
me60042ad <- rowMeans(e60042ad)
e60042adh <- e60042ad[me60042ad>=lowexpthr,]

# next, we reduce genes by variance
# to avoid lowly varied genes from correlation based analysis
lowvarthr = 0.5
de60042adhn <- apply(e60042adh[,e60042adns],1,sd)
de60042adht <- apply(e60042adh[,e60042adts],1,sd)
de60042adh <- apply(cbind(de60042adhn,de60042adht),1,max)
e60042adhv <- e60042adh[de60042adh>=lowvarthr,]
print(paste(nrow(e60042adhv),' genes are selected for GSe60042 adenoma'))

remove(e60042ad)
remove(e60042adh)

# finally, calculate spearman correlations and 
# multiple testing corrected p-values
library(Hmisc)
library(metaSEM)

# Calculation of degree values for the cancer networks
# self correlations are set to 0
# positive and negative correlations are considered separately
# connection requires significant correlation based on randomization
lowthr = 0.2
highthr = 0.5

# real correlations based on real normal and tumor sample sets
# for e8671 (adenoma), samples (columns) 1-32  are normal, 33-64 are tumor
# for e18105 (carcinoma 1), samples (columns) 1-17 are normal, 18-34 are tumor
noradeexp <- e60042adhv[,e60042adns];tumadeexp <- e60042adhv[,e60042adts]

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

# randomized samples based correlations
# take random sample sets of same size
# skip first columns which is IDs
for (jj in 1:norands){
  print(jj)
  
  sc = c(1:e60042adts[length(e60042adts)])
  randset = sample(sc,e60042adns[length(e60042adns)])
  randset2 = sc[-randset]
  randnoradeexp <- e60042adhv[,randset];randtumadeexp <- e60042adhv[,randset2]
  
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

# raposn, rapost, .... capture random connections so they will be substracted  
# posn2, post2, negn2 and negt2 ignore the correlation on the other sample type
# differential correlation: hardposn, hardpost, hardnegn and hardnegt
# correlation >= highthr in the sample and lower than 0.05 in the randomized correlations
# and <= lowthr in corresponding normal or tumor

aposn2 = aposn2 - raposn
hardposnx = aposn2 - apost1
hardposnx[hardposnx < 0] = 0

apost2 = apost2 - rapost
hardpostx = apost2 - aposn1
hardpostx[hardpostx < 0] = 0

anegn2 = anegn2 - ranegn
hardnegnx = anegn2 - anegt1
hardnegnx[hardnegnx < 0] = 0

anegt2 = anegt2 - ranegt
hardnegtx = anegt2 - anegn1
hardnegtx[hardnegtx < 0] = 0

napd <- rowSums(hardposnx);nand <- rowSums(hardnegnx)
tapd <- rowSums(hardpostx);tand <- rowSums(hardnegtx)

# less and more connected genes based on their weighted degree difference
difapd <- tapd-napd
lessapd <- difapd[difapd<0];lessapd <- abs(lessapd)
moreapd <- difapd[difapd>0]


difand <- tand-nand
lessand <- difand[difand<0];lessand <- abs(lessand)
moreand <- difand[difand>0]


# write less and more connected genes with weights based on weighted degree difference
write.table(lessapd,file='lessapd_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(moreapd,file='moreapd_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessand,file='lessand_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(moreand,file='moreand_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)


# write top ranking genes
lessapd5 = lessapd[lessapd >= quantile(lessapd,probs=0.95)]
lessapd10 = lessapd[lessapd >= quantile(lessapd,probs=0.9)]
lessapd20 = lessapd[lessapd>=quantile(lessapd,probs=0.8)] 
write.table(lessapd5,file='lessapd_top5percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessapd10,file='lessapd_top10percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessapd20,file='lessapd_top20percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)

perc5 = quantile(moreapd,probs=0.95)
perc10 = quantile(moreapd,probs=0.9)
perc20 = quantile(moreapd,probs=0.8)
write.table(moreapd[moreapd>=perc5],file='moreapd_top5percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(moreapd[moreapd>=perc10],file='moreapd_top10percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(moreapd[moreapd>=perc20],file='moreapd_top20percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)

perc5 = quantile(lessand,probs=0.95)
perc10 = quantile(lessand,probs=0.9)
perc20 = quantile(lessand,probs=0.8)
write.table(lessand[lessand>=perc5],file='lessand_top5percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessand[lessand>=perc10],file='lessand_top10percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessand[lessand>=perc20],file='lessand_top20percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)

perc5 = quantile(moreand,probs=0.95)
perc10 = quantile(moreand,probs=0.9)
perc20 = quantile(moreand,probs=0.8)
write.table(moreand[moreand>=perc5],file='moreand_top5percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(moreand[moreand>=perc10],file='moreand_top10percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(moreand[moreand>=perc20],file='moreand_top20percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)


# neighbors of top ranking genes
# write top ranking genes
perc5 = quantile(lessapd,probs=0.95)
perc10 = quantile(lessapd,probs=0.9)
perc20 = quantile(lessapd,probs=0.8)
write.table(lessapd[lessapd>=perc5],file='lessapd_top5percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessapd[lessapd>=perc10],file='lessapd_top10percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)
write.table(lessapd[lessapd>=perc20],file='lessapd_top20percent_Adrenal.txt',row.names=TRUE,col.names=FALSE,sep=',',quote=FALSE)



# plot degs and deg differences
png(file='adeno_carci_degs_Adrenal.png',width=240,height=600)
par(mfrow = c(2,1))

# plot ABSOLUTE degree values
boxplot(napd,tapd,main='positive adenoma',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5), xaxt='n',
        cex.axis=1.5,cex.main=1.75,cex.lab=1.5)
axis(1,at=c(0.6,1.5),labels=c('normal','tumor'),tick=FALSE,cex.axis=1.5)


boxplot(nand,tand,main='negative adenoma',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5), xaxt='n',
        cex.axis=1.5,cex.main=1.75,cex.lab=1.5)
axis(1,at=c(0.6,1.5),labels=c('normal','tumor'),tick=FALSE,cex.axis=1.5)

dev.off()

# mann whitney tests for comparing degree dist. of tumor to normal (paired)
wt <- wilcox.test(napd,tapd,paired=TRUE,alternative='greater')
print(paste('adenoma positive paired wilcox test p-value:',wt$p.value))
wt <- wilcox.test(nand,tand,paired=TRUE,alternative='greater')
print(paste('adenoma negative paired wilcox test p-value:',wt$p.value))

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


