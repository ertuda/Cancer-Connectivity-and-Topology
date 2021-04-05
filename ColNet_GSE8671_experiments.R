# GSE8671 experiments
# Experimenting with various options, cutoffs ...
# By Ertugrul Dalgic, PhD
# 2015-2018

# index of normal and tumor samples
e8671ns = 1:32
e8671ts = 33:64

# Normalize the raw dataset with RMA
library(affy)
setwd('C:/Users/Dell/Documents/Projects/Colorectal Topology/GSE8671/')
eade <- read.AnnotatedDataFrame('GSE8671_targets.txt',header=FALSE,
                                row.names=1,as.is=TRUE)
raweade = ReadAffy(filenames=pData(eade)$FileName,phenoData=eade)
normeade = rma(raweade)
e8671 = exprs(normeade)
remove(normeade);remove(raweade)

# first reducing affy ids
# collapse by max var
library(hgu133plus2.db)
mapped_probes <- mappedkeys(hgu133plus2ENTREZID)
hid <- as.list(hgu133plus2ENTREZID[mapped_probes])
e8671gid <- e8671[intersect(rownames(e8671),names(hid)),]
rownames(e8671gid) <- sapply(rownames(e8671gid),function(i) hid[[i]])
# use var to select the probe set
ax <- apply(e8671gid,1,sd)
sax <- split(ax,as.factor(names(ax)))
ssax <- sapply(sax,max)
sind <- sapply(seq_along(ax),function(i) ax[i]==ssax[names(ax)[i]])
e8671j <- e8671gid[sind,]
remove(e8671gid)
remove(e8671)

# next, we reduce genes by average expression level
# to avoid lowly expressed genes from correlation based analysis
lowexpthr = 5
me8671j <- rowMeans(e8671j)
#png('GSE8671_avgexp.png')
#par(mfrow=c(1,1))
#hist(me8671j,breaks=40,xlim=c(0,15),ylim=c(0,2500),xlab='Average expression',
#     cex.axis=1.25,cex.lab=1.25)
#dev.off()
e8671jh <- e8671j[me8671j>=lowexpthr,]

# next, we reduce genes by quartile coefficient of dispersion (~variance)
# to avoid lowly varied genes from correlation based analysis
# get genes with a variance above threshold in either normal or tumor samples
lowvarthr = 0.035
de8671jhn <- apply(e8671jh[,1:32],1,function(i) ((quantile(i,0.75)-quantile(i,0.25))/
                                             (quantile(i,0.75)+quantile(i,0.25))))
de8671jht <- apply(e8671jh[,33:64],1,function(i) ((quantile(i,0.75)-quantile(i,0.25))/
                                                  (quantile(i,0.75)+quantile(i,0.25))))

#png('GSE8671_var.png')
#par(mfrow=c(1,1))
#hist(de8671jhn,breaks=40,xlim=c(0,0.4),xlab="Coefficient
#     of variation",cex.axis=1.25,cex.lab=1.25)
#dev.off()
e8671jhvn <- e8671jh[de8671jhn>=lowvarthr,]
e8671jhvt <- e8671jh[de8671jht>=lowvarthr,]
e8671jhv <- e8671jh[intersect(rownames(e8671jhvn),rownames(e8671jhvt)),]


# next, we reduce genes by variance
# to avoid lowly varied genes from correlation based analysis
lowvarthr = 0.5
de8671jhn <- apply(e8671jh[,e8671ns],1,sd)
de8671jht <- apply(e8671jh[,e8671ts],1,sd)
de8671jh <- apply(cbind(de8671jhn,de8671jht),1,max)
e8671jhv <- e8671jh[de8671jh>=lowvarthr,]


print(paste(nrow(e8671jhv),' genes are selected for GSE8671'))

# finally, calculate spearman correlations and 
# multiple testing corrected p-values
library(Hmisc)
library(metaSEM)

# finally, calculate correlations based on the real and random gene sets 
# Calculation of weighted degree values for the cancer networks
# self correlations are set to 0
# positive and negative correlations are considered separately
# connection requires correlation 
lowthr = 0.2
highthr = 0.8

nexp <- e8671jhv[,e8671ns]
ncor <- rcorr(t(nexp),type='pearson')

posn1 <- ncor$r
diag(posn1) = 0
posn1[posn1 < lowthr] = 0
posn1[posn1 != 0] = 1
negn1 <- ncor$r
diag(negn1) = 0
negn1[negn1 > -lowthr] = 0
negn1[negn1 != 0] = 1

posn2 <- ncor$r
diag(posn2) = 0
posn2[posn2 < highthr] = 0
posn2[posn2 != 0] = 1
negn2 <- ncor$r
diag(negn2) = 0
negn2[negn2 > -highthr] = 0
negn2[negn2 != 0] = 1

texp <- e8671jhv[,e8671ts]
tcor <- rcorr(t(texp),type='pearson')

post1 <- tcor$r
diag(post1) = 0
post1[post1 < lowthr] = 0
post1[post1 != 0] = 1
negt1 <- tcor$r
diag(negt1) = 0
negt1[negt1 > -lowthr] = 0
negt1[negt1 != 0] = 1

post2 <- tcor$r
diag(post2) = 0
post2[post2 < highthr] = 0
post2[post2 != 0] = 1
negt2 <- tcor$r
diag(negt2) = 0
negt2[negt2 > -highthr] = 0
negt2[negt2 != 0] = 1

# posn2, post2, negn2 and negt2 ignore the correlation on the other sample type
# differential correlation: hardposn, hardpost, hardnegn and hardnegt
# correlation >= highthr in the sample and <= lowthr in corresponding normal or tumor

hardposn = posn2 - post1
hardposn[hardposn == -1] = 0

hardpost = post2 - posn1
hardpost[hardpost == -1] = 0

hardnegn = negn2 - negt1
hardnegn[hardnegn == -1] = 0

hardnegt = negt2 - negn1
hardnegt[hardnegt == -1] = 0

# degree values as sum of weighted degree values
# normal
pnd <- rowSums(posn2);nnd <- rowSums(negn2)
ptd <- rowSums(post2);ntd <- rowSums(negt2)

# differential
pnd <- rowSums(hardposn);nnd <- rowSums(hardnegn)
ptd <- rowSums(hardpost);ntd <- rowSums(hardnegt)

# mann whitney tests for comparing degree dist. of tumor to normal (paired)
wt <- wilcox.test(pnd,ptd,paired=TRUE,alternative='greater')
print(paste('adenoma positive paired wilcox test p-value:',wt$p.value))
wt <- wilcox.test(nnd,ntd,paired=TRUE,alternative='greater')
print(paste('adenoma negative paired wilcox test p-value:',wt$p.value))


# plot degs and deg differences
tiff(file='adeno_degs_exp_cor08.png.tiff',width=3.3, height=3,units="in",res=600)
par(mfrow = c(1,2),mar=c(2.3,3.9,2.3,1))

# plot ABSOLUTE degree values
boxplot(pnd,ptd,main='positive',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5), xaxt='n',
        cex.axis=0.7,cex.main=1,cex.lab=1)
axis(1,at=c(0.6,1.5),labels=c('N','T'),tick=FALSE,cex.axis=1)
boxplot(nnd,ntd,main='negative',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5), xaxt='n',
        cex.axis=0.7,cex.main=1,cex.lab=1)
axis(1,at=c(0.6,1.5),labels=c('N','T'),tick=FALSE,cex.axis=1)
dev.off()




# plot ABSOLUTE degree values
boxplot(napd,tapd,main='Ade 2 Pos',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5), xaxt='n',
        cex.axis=1.15,cex.lab=1.25)
axis(1,at=c(0.6,1.5),labels=c('N','T'),tick=FALSE,cex.axis=1.15)

boxplot(ncpd,tcpd,main='Car 2 Pos',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5), xaxt='n',
        cex.axis=1.15,cex.lab=1.25)
axis(1,at=c(0.6,1.5),labels=c('N','T'),tick=FALSE,cex.axis=1.15)

boxplot(nand,tand,main='Ade 2 Neg',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5), xaxt='n',
        cex.axis=1.15,cex.lab=1.25)
axis(1,at=c(0.6,1.5),labels=c('N','T'),tick=FALSE,cex.axis=1.15)

boxplot(ncnd,tcnd,main='Car 2 Neg',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5), xaxt='n',
        cex.axis=1.15,cex.lab=1.25)
axis(1,at=c(0.6,1.5),labels=c('N','T'),tick=FALSE,cex.axis=1.15)
dev.off()




###################
# plot smaller sets
png(file='adeno_degs_smallersets.png',width=480,height=1000)
par(mfrow = c(2,1))

# plot ABSOLUTE degree values
boxplot(pnd5,ptd5,pnd10,ptd10,pnd15,ptd15,pnd20,ptd20,pnd25,ptd25,pnd30,ptd30,pnd,ptd,
        main='positive',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5,2.5,3.4,4.4,5.3,6.3,7.2,8.2,9.1,10,10.9,11.8,12.7), 
        xaxt='n', cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
legend('topright',legend=c('normal','tumor'),fill=c('white','gray50'))
axis(1,at=c(1.1,3,4.9,6.8,8.7,10.5,12.4),labels=c('5','10','15','20','25','30','32'),
     tick=FALSE,cex.axis=1.25)
boxplot(nnd5,ntd5,nnd10,ntd10,nnd15,ntd15,nnd20,ntd20,nnd25,ntd25,nnd30,ntd30,nnd,ntd,
        main='negative',ylab='degree values',
        col=c('white','gray50'), at=c(0.6,1.5,2.5,3.4,4.4,5.3,6.3,7.2,8.2,9.1,10,10.9,11.8,12.7), 
        xaxt='n', cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
legend('topright',legend=c('normal','tumor'),fill=c('white','gray50'))
axis(1,at=c(1.1,3,4.9,6.8,8.7,10.5,12.4),labels=c('5','10','15','20','25','30','32'),
     tick=FALSE,cex.axis=1.25)
dev.off()



##################
# continue ...
# randomization !! 
############################################
# take random sample sets of size 17
for (jj in 1:100){
    print(jj)
    normalset <- sample(1:32,17)
    tumorset <- normalset + 32
    norcarexp <- e8671jhv[,normalset]; tumcarexp <- e8671jhv[,tumorset]

# # test the original sample sets    
# for (jj in 1:1){
#     norcarexp <- e8671jhv[,1:32]; tumcarexp <- e8671jhv[,33:64]
    
    # SPEARMAN BASED CORRELATIONS
    # normal carcinoma
    ncorcar <- rcorr(t(norcarexp),type='spearman')
    ncorvalcar <- sm2vec(ncorcar$r)
    ncorsigcar <- sm2vec(ncorcar$P)
    ncorsigmccar <- p.adjust(ncorsigcar,method='BH')
    ncorsigmcmatcar <- vec2sm(ncorsigmccar)
    # tumor carcinoma
    tcorcar <- rcorr(t(tumcarexp),type='spearman')
    tcorvalcar <- sm2vec(tcorcar$r)
    tcorsigcar <- sm2vec(tcorcar$P)
    tcorsigmccar <- p.adjust(tcorsigcar,method='BH')
    tcorsigmcmatcar <- vec2sm(tcorsigmccar)
    
    # connection requires spearman correlations p value <= pvthr
    norcarpos <- ncorcar$r
    diag(norcarpos) = 0;norcarpos[norcarpos<0] = 0;norcarpos[ncorsigmcmatcar>pvthr] = 0
    norcarneg <- ncorcar$r
    diag(norcarneg) = 0;norcarneg[norcarneg>0] = 0;norcarneg[ncorsigmcmatcar>pvthr] = 0
    norcarneg <- abs(norcarneg)
    tumcarpos <- tcorcar$r
    diag(tumcarpos) = 0;tumcarpos[tumcarpos<0] = 0;tumcarpos[tcorsigmcmatcar>pvthr] = 0
    tumcarneg <- tcorcar$r
    diag(tumcarneg) = 0;tumcarneg[tumcarneg>0] = 0;tumcarneg[tcorsigmcmatcar>pvthr] = 0
    tumcarneg <- abs(tumcarneg)
    
#     # PEARSON BASED CORRELATIONS
#     # normal carcinoma
#     ncorcarpr <- rcorr(t(norcarexp),type='pearson')
#     ncorvalcarpr <- sm2vec(ncorcarpr$r)
#     ncorsigcarpr <- sm2vec(ncorcarpr$P)
#     ncorsigmccarpr <- p.adjust(ncorsigcarpr,method='BH')
#     ncorsigmcmatcarpr <- vec2sm(ncorsigmccarpr)
#     # tumor carcinoma
#     tcorcarpr <- rcorr(t(tumcarexp),type='pearson')
#     tcorvalcarpr <- sm2vec(tcorcarpr$r)
#     tcorsigcarpr <- sm2vec(tcorcarpr$P)
#     tcorsigmccarpr <- p.adjust(tcorsigcarpr,method='BH')
#     tcorsigmcmatcarpr <- vec2sm(tcorsigmccarpr)
#     
#     # connection requires pearson correlations p value <= pvthr
#     norcarpos <- ncorcarpr$r
#     diag(norcarpos) = 0;norcarpos[norcarpos<0] = 0;norcarpos[ncorsigmcmatcarpr>pvthr] = 0
#     norcarneg <- ncorcarpr$r
#     diag(norcarneg) = 0;norcarneg[norcarneg>0] = 0;norcarneg[ncorsigmcmatcarpr>pvthr] = 0
#     norcarneg <- abs(norcarneg)
#     tumcarpos <- tcorcarpr$r
#     diag(tumcarpos) = 0;tumcarpos[tumcarpos<0] = 0;tumcarpos[tcorsigmcmatcarpr>pvthr] = 0
#     tumcarneg <- tcorcarpr$r
#     diag(tumcarneg) = 0;tumcarneg[tumcarneg>0] = 0;tumcarneg[tcorsigmcmatcarpr>pvthr] = 0
#     tumcarneg <- abs(tumcarneg)
    
    # degree values 
    ncpd <- rowSums(norcarpos);ncnd <- rowSums(norcarneg)
    tcpd <- rowSums(tumcarpos);tcnd <- rowSums(tumcarneg)
    
    norposdegmeds2 <- c(norposdegmeds2,quantile(ncpd,0.5))
    norposdegfqs2 <- c(norposdegfqs2,quantile(ncpd,0.25))
    norposdegtqs2 <- c(norposdegtqs2,quantile(ncpd,0.75))
    tumposdegmeds2 <- c(tumposdegmeds2,quantile(tcpd,0.5))
    tumposdegfqs2 <- c(tumposdegfqs2,quantile(tcpd,0.25))
    tumposdegtqs2 <- c(tumposdegtqs2,quantile(tcpd,0.75))
    nornegdegmeds2 <- c(nornegdegmeds2,quantile(ncnd,0.5))
    nornegdegfqs2 <- c(nornegdegfqs2,quantile(ncnd,0.25))
    nornegdegtqs2 <- c(nornegdegtqs2,quantile(ncnd,0.75))
    tumnegdegmeds2 <- c(tumnegdegmeds2,quantile(tcnd,0.5))
    tumnegdegfqs2 <- c(tumnegdegfqs2,quantile(tcnd,0.25))
    tumnegdegtqs2 <- c(tumnegdegtqs2,quantile(tcnd,0.75))
}

# plot degree medians along the sample size range from 5 to 32
png(file='GSE8671_samplesizes_effect.png',width = 760, height = 480)
par(mfrow=c(1,2))
maxtqs <- c(max(norposdegtqs2),max(tumposdegtqs2)); maxboth <- max(maxtqs)
plot(norposdegmeds2,type='l',col='blue',lwd=4,ylim=c(0,maxboth+5),main='positive',
     xlim=c(5,32+5),ylab='degree values',xlab='sample size',cex=1.5)
lines(norposdegfqs2,lwd=0.5,lty=3,col='lightblue')
lines(norposdegtqs2,lwd=0.5,lty=3,col='lightblue')
lines(tumposdegmeds2,lwd=4,col='red')
lines(tumposdegfqs2,lwd=0.5,lty=3,col='pink')
lines(tumposdegtqs2,lwd=0.5,lty=3,col='pink')
legend('topright',legend=c('normal','tumor'),fill=c('blue','red'))

maxtqs <- c(max(nornegdegtqs2),max(tumnegdegtqs2)); maxboth <- max(maxtqs)
plot(nornegdegmeds2,type='l',col='blue',lwd=4,ylim=c(0,maxboth+5),main='negative',
     xlim=c(5,32+5),ylab='degree values',xlab='sample size',cex=1.5)
lines(nornegdegfqs2,lwd=0.5,lty=3,col='lightblue')
lines(nornegdegtqs2,lwd=0.5,lty=3,col='lightblue')
lines(tumnegdegmeds2,lwd=4,col='red')
lines(tumnegdegfqs2,lwd=0.5,lty=3,col='pink')
lines(tumnegdegtqs2,lwd=0.5,lty=3,col='pink')
legend('topright',legend=c('normal','tumor'),fill=c('blue','red'))
dev.off()

# plot difference of degree medians for 100 randomizations
png(file='GSE8671_randomsamplesets_effect.png')
par(mfrow=c(1,1))
hist((tumposdegmeds2-norposdegmeds2),(tumnegdegmeds2-nornegdegmeds2)
     ,ylab='difference of tumor from normal',xlab='sample size',cex=1.5)
axis(1,at=c(1,2),labels=c('positive','negative'))
dev.off()

png(file='GSE8671_degs_jetset_maxvar.png',width = 760, height = 480)
par(mfrow=c(1,2))
boxplot(ncpdj,tcpdj,ncpdm,tcpdm,main='positive',ylab='degree values',
        col=c('lightcyan3','firebrick3'), at=c(0.6,1.5,2.9,3.8), xaxt='n',
        cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
legend('topright',legend=c('normal','tumor'),fill=c('lightcyan3','firebrick3'))
axis(1,at=c(1,3.4),labels=c('Jetset','Max. var.'),tick=FALSE,cex.axis=1.25)
boxplot(ncndj,tcndj,ncndm,tcndm,main='negative',ylab='degree values',
        col=c('lightcyan3','firebrick3'), at=c(0.6,1.5,2.9,3.8), xaxt='n',
        cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
legend('topright',legend=c('normal','tumor'),fill=c('lightcyan3','firebrick3'))
axis(1,at=c(1,3.4),labels=c('Jetset','Max. var.'),tick=FALSE,cex.axis=1.25)
dev.off()

png(file='GSE8671_degs_expthrs_4_5_6.png',width = 760, height = 480)
par(mfrow=c(1,2))
boxplot(ncpd4,tcpd4,ncpd5,tcpd5,ncpd6,tcpd6,main='positive',ylab='degree values',
        col=c('lightcyan3','firebrick3'), at=c(0.6,1.5,2.9,3.8,5.2,6.1), xaxt='n',
        cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
legend('topright',legend=c('normal','tumor'),fill=c('lightcyan3','firebrick3'))
axis(1,at=c(1,3.4,5.7),labels=c('4','5','6'),tick=FALSE,cex.axis=1.25)
boxplot(ncnd4,tcnd4,ncnd5,tcnd5,ncnd6,tcnd6,main='negative',ylab='degree values',
        col=c('lightcyan3','firebrick3'), at=c(0.6,1.5,2.9,3.8,5.2,6.1), xaxt='n',
        cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
legend('topright',legend=c('normal','tumor'),fill=c('lightcyan3','firebrick3'))
axis(1,at=c(1,3.4,5.7),labels=c('4','5','6'),tick=FALSE,cex.axis=1.25)
dev.off()

png(file='GSE8671_degs_varthrs_04_05_06.png',width = 760, height = 480)
par(mfrow=c(1,2))
boxplot(ncpd04,tcpd04,ncpd05,tcpd05,ncpd06,tcpd06,main='positive',ylab='degree values',
        col=c('lightcyan3','firebrick3'), at=c(0.6,1.5,2.9,3.8,5.2,6.1), xaxt='n',
        cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
legend('topright',legend=c('normal','tumor'),fill=c('lightcyan3','firebrick3'))
axis(1,at=c(1,3.4,5.7),labels=c('0.04','0.05','0.06'),tick=FALSE,cex.axis=1.25)
boxplot(ncnd04,tcnd04,ncnd05,tcnd05,ncnd06,tcnd06,main='negative',ylab='degree values',
        col=c('lightcyan3','firebrick3'), at=c(0.6,1.5,2.9,3.8,5.2,6.1), xaxt='n',
        cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
legend('topright',legend=c('normal','tumor'),fill=c('lightcyan3','firebrick3'))
axis(1,at=c(1,3.4,5.7),labels=c('0.04','0.05','0.06'),tick=FALSE,cex.axis=1.25)
dev.off()

png(file='GSE8671_degs_spearman_pearson.png',width = 760, height = 480)
par(mfrow=c(1,2))
boxplot(ncpds,tcpds,ncpdp,tcpdp,main='positive',ylab='degree values',
        col=c('lightcyan3','firebrick3'), at=c(0.6,1.5,2.9,3.8), xaxt='n',
        cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
legend('topright',legend=c('normal','tumor'),fill=c('lightcyan3','firebrick3'))
axis(1,at=c(1,3.4),labels=c('Spearman','Pearson'),tick=FALSE,cex.axis=1.25)
boxplot(ncnds,tcnds,ncndp,tcndp,main='negative',ylab='degree values',
        col=c('lightcyan3','firebrick3'), at=c(0.6,1.5,2.9,3.8), xaxt='n',
        cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
legend('topright',legend=c('normal','tumor'),fill=c('lightcyan3','firebrick3'))
axis(1,at=c(1,3.4),labels=c('Spearman','Pearson'),tick=FALSE,cex.axis=1.25)
dev.off()

png(file='GSE8671_degs_Spearmanthrs_001_005_01.png',width = 760, height = 480)
par(mfrow=c(1,2))
boxplot(ncpd01,tcpd01,ncpd05,tcpd05,ncpd1,tcpd1,main='positive',ylab='degree values',
        col=c('lightcyan3','firebrick3'), at=c(0.6,1.5,2.9,3.8,5.2,6.1), xaxt='n',
        cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
legend('topright',legend=c('normal','tumor'),fill=c('lightcyan3','firebrick3'))
axis(1,at=c(1,3.4,5.7),labels=c('0.01','0.05','0.1'),tick=FALSE,cex.axis=1.25)
boxplot(ncnd01,tcnd01,ncnd05,tcnd05,ncnd1,tcnd1,main='negative',ylab='degree values',
        col=c('lightcyan3','firebrick3'), at=c(0.6,1.5,2.9,3.8,5.2,6.1), xaxt='n',
        cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
legend('topright',legend=c('normal','tumor'),fill=c('lightcyan3','firebrick3'))
axis(1,at=c(1,3.4,5.7),labels=c('0.01','0.05','0.1'),tick=FALSE,cex.axis=1.25)
dev.off()


# plot degs and deg differences
png(file='GSE8671_degs.png')
par(mfrow = c(2,2))
# plot ABSOLUTE degree values
boxplot(ncpd,tcpd,main='positive',ylab='degree values',
        col=c('lightcyan3','firebrick3'),cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
axis(1,at=c(1,2),labels=c('normal','tumor'),tick=FALSE,cex.axis=1.25)
boxplot(ncnd,tcnd,main='negative',ylab='degree values',
        col=c('lightcyan3','firebrick3'),cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
axis(1,at=c(1,2),labels=c('normal','tumor'),tick=FALSE,cex.axis=1.25)
# plot DIFFERENCE of degree values
hist(tcpd-ncpd,main='carcinoma positive',col='magenta',xlab='difference of degree',
     ylab='frequency',cex.axis=1.25,cex.main=1.25,cex.lab=1.25,xlim=c(-90,90))
abline(v=0,lty=2,lwd=4)
hist(tcnd-ncnd,main='carcinoma negative',col='indianred4',xlab='difference of degree',
     ylab='frequency',cex.axis=1.25,cex.main=1.25,cex.lab=1.25,xlim=c(-40,40))
abline(v=0,lty=2,lwd=4)
dev.off()

######################

# DIFFERENTIAL NETWORKS
# in order to construct differential network difference of cor. p-values
# are considered; the threshold cutting above can create unwanted list of 
# pairs (a little above the thr and and a little down in the paired samples)
# also a significant p val for both spearman and pearson is required 
# to get very significant pairs

# use pearson based cor. p values as a condition in the following step
ncsppr <- ncorsigmcmatcarpr; diag(ncsppr)=10
ncsppr[ncorcarpr$r<0]=10; ncsppr[ncsppr>pvthr]=10

ncsnpr <- ncorsigmcmatcarpr; diag(ncsnpr)=10
ncsnpr[ncorcarpr$r>0]=10; ncsnpr[ncsnpr>pvthr]=10

tcsppr <- tcorsigmcmatcarpr; diag(tcsppr)=10
tcsppr[tcorcarpr$r<0]=10; tcsppr[tcsppr>pvthr]=10

tcsnpr <- tcorsigmcmatcarpr; diag(tcsnpr)=10
tcsnpr[tcorcarpr$r>0]=10; tcsnpr[tcsnpr>pvthr]=10

ncsp2pr <- ncorsigmcmatcarpr; diag(ncsppr)=10;ncsppr[ncorcarpr$r<0]=10
ncsn2pr <- ncorsigmcmatcarpr; diag(ncsnpr)=10;ncsnpr[ncorcarpr$r>0]=10
tcsp2pr <- tcorsigmcmatcarpr; diag(tcsppr)=10;tcsppr[tcorcarpr$r<0]=10
tcsn2pr <- tcorsigmcmatcarpr; diag(tcsnpr)=10;tcsnpr[tcorcarpr$r>0]=10

# use spearman based cor. p values to get the pairs
ncsp <- ncorsigmcmatcar; diag(ncsp)=10;ncsp[ncorcar$r<0]=10
ncsp[ncsp>pvthr]=10; ncsp[ncsppr>pvthr]=10

ncsn <- ncorsigmcmatcar; diag(ncsn)=10;ncsn[ncorcar$r>0]=10
ncsn[ncsn>pvthr]=10; ncsn[ncsnpr>pvthr]=10

tcsp <- tcorsigmcmatcar; diag(tcsp)=10;tcsp[tcorcar$r<0]=10
tcsp[tcsp>pvthr]=10; tcsp[tcsppr>pvthr]=10

tcsn <- tcorsigmcmatcar; diag(tcsn)=10;tcsn[tcorcar$r>0]=10
tcsn[tcsn>pvthr]=10; tcsn[tcsnpr>pvthr]=10

ncsp2 <- ncorsigmcmatcar; diag(ncsp)=10;ncsp[ncorcar$r<0]=10
ncsn2 <- ncorsigmcmatcar; diag(ncsn)=10;ncsn[ncorcar$r>0]=10
tcsp2 <- tcorsigmcmatcar; diag(tcsp)=10;tcsp[tcorcar$r<0]=10
tcsn2 <- tcorsigmcmatcar; diag(tcsn)=10;tcsn[tcorcar$r>0]=10

# get the min of paired p values in order to avoid any possible correlation
# in the paired sample set 
ncspmin <- pmin(ncsp2,ncsp2pr); ncsnmin <- pmin(ncsn2,ncsn2pr)
tcspmin <- pmin(tcsp2,tcsp2pr); tcsnmin <- pmin(tcsn2,tcsn2pr)

ncspdif <- tcspmin-ncsp; ncsndif <- tcsnmin-ncsn
tcspdif <- ncspmin-tcsp; tcsndif <- ncsnmin-tcsn

sigposnc <- sm2vec(ncspdif); sigposnc <- intersect(sigposnc[sigposnc>0],sigposnc[sigposnc<1])
sigpostc <- sm2vec(tcspdif); sigpostc <- intersect(sigpostc[sigpostc>0],sigpostc[sigpostc<1])

signegnc <- sm2vec(ncsndif); signegnc <- intersect(signegnc[signegnc>0],signegnc[signegnc<1])
signegtc <- sm2vec(tcsndif); signegtc <- intersect(signegtc[signegtc>0],signegtc[signegtc<1])

# plot distribution of pairs with a difference of p-value from the paired normal/tumor
png(file='GSE44076_adecar_pvaldiff.png',width = 800, height = 800)
par(mfrow = c(2,1))
plot(density(sigpostc),col='red',lwd=3,main='positive carcinoma',
     xlab='p val. diff.',cex.axis=1.3,cex.main=1.3,cex.lab=1.3,xlim=c(-0.5,1.2))
lines(density(sigposnc),col='blue',lwd=3)
legend('topright',legend=c('tumor','normal'),fill=c('red','blue'))
plot(density(signegnc),col='blue',lwd=3,main='negative carcinoma',
     xlab='p val. diff.',cex.axis=1.3,cex.main=1.3,cex.lab=1.3,ylim=c(0,5))
lines(density(signegtc),col='red',lwd=3)
legend('topright',legend=c('tumor','normal'),fill=c('red','blue'))
dev.off()




#########################
# write carcinoma gene pairs with their connection weights and export as text file
expfname = 'GSE44076_adenoma_carcinoma_pospairs_verysig.txt'
pvaldifthr = 0.6
carpospairs <- c('Gene1','Gene2','Weight')
genesnoname <- setdiff(rownames(ncorcar$P),aidgid$EntrezID)
for (i in 1:(nrow(ncspdif)-1)){
    for (j in (i+1):nrow(ncspdif)){
        if (ncspdif[i,j] > pvaldifthr & ncspdif[i,j] < 1){
            if (any(genesnoname == rownames(ncorcar$P)[i]) == TRUE){
                namei = rownames(ncorcar$P)[i]}
            else{
                namei <- aidgid[aidgid$EntrezID==rownames(ncorcar$P)[i],'symbol']}
            if (any(genesnoname == rownames(ncorcar$P)[i]) == TRUE){
                namej = rownames(ncorcar$P)[j]}
            else{
                namej <- aidgid[aidgid$EntrezID==rownames(ncorcar$P)[j],'symbol']}
            carpospairs <- rbind(carpospairs,c(namei,namej,abs(ncspdif[i,j])))
        }}}
write.table(carpospairs,file=expfname,row.names=FALSE,col.names=FALSE,sep='\t')

# plot the carcinoma expression values of a gene pair of interest 
givengene1 = 'SGK1'
givengene2 = 'MAOA'
par(mfrow=c(1,2))
plot(norcarexp[aidgid[aidgid$symbol==givengene1,'EntrezID'],],
     norcarexp[aidgid[aidgid$symbol==givengene2,'EntrezID'],],
     xlab=givengene1,ylab=givengene2,main='adenoma',cex.lab=1.3,cex.main=1.3)
plot(tumcarexp[aidgid[aidgid$symbol==givengene1,'EntrezID'],],
     tumcarexp[aidgid[aidgid$symbol==givengene2,'EntrezID'],],
     xlab=givengene1,ylab=givengene2,main='carcinoma',cex.lab=1.3,cex.main=1.3)
print('adenoma spearman correlation test')
cor.test(norcarexp[aidgid[aidgid$symbol==givengene1,'EntrezID'],],
         norcarexp[aidgid[aidgid$symbol==givengene2,'EntrezID'],],method='spearman')
print('carcinoma spearman correlation test')
cor.test(tumcarexp[aidgid[aidgid$symbol==givengene1,'EntrezID'],],
         tumcarexp[aidgid[aidgid$symbol==givengene2,'EntrezID'],],method='spearman')



