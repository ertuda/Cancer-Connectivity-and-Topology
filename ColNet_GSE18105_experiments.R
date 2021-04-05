# GSE18105 experiments
# Experimenting with various options, cutoffs ...
# By Ertugrul Dalgic, PhD
# 2015-2018

# Normalize the raw datasets with RMA from affy library
library(affy)

# Colorectal Carcinoma Dataset (GSE18105)
setwd('C:/Users/Dell/Documents/Projects/Colorectal Topology/GSE18105/')
ecar <- read.AnnotatedDataFrame('GSE18105_targets.txt',header=FALSE,
                                row.names=1,as.is=TRUE)
rawecar = ReadAffy(filenames=pData(ecar)$FileName,phenoData=ecar)
normecar = rma(rawecar)
e18105 = exprs(normecar)
remove(normecar);remove(rawecar)

setwd('C:/Users/Dell/Documents/Projects/Colorectal Topology/')

# first reducing affy ids based on the jetset package
# match affy ids to gene ids; get the BEST probesets (best='TRUE')
# name-id matching is in aidgid; it can be used later in the network
library(jetset)
jetsetfile = 'jetset.scores.hgu133plus2_3.4.0.csv'
aidgid <- read.csv(jetsetfile,colClasses='character')
aidgid <- aidgid[aidgid$best=='TRUE',]
print(paste(length(intersect(aidgid$probeset,rownames(e18105))),'of ',
            length(aidgid$probeset),'affy ids are present in the expression datasets'))
ind18105 <- as.vector(sapply(aidgid$EntrezID,
                             function(i) which(rownames(e18105)==jmap('hgu133plus2',eg=i))))
e18105j <- e18105[ind18105,]
rownames(e18105j) <- sapply(rownames(e18105j),
                            function(i) aidgid[which(aidgid$probeset==i),'EntrezID'])

remove(e18105)

# next, we reduce genes by average expression level
# to avoid lowly expressed genes from correlation based analysis
lowexpthr = 5
me18105j <- rowMeans(e18105j)
# png('GSE18105_avgexp.png')
# par(mfrow=c(1,1))
# hist(me18105j,breaks=40,xlim=c(0,15),ylim=c(0,2500),xlab='Average expression',
#      cex.axis=1.25,cex.lab=1.25)
# dev.off()
e18105jh <- e18105j[me18105j>=lowexpthr,]

# next, we reduce genes by coefficient of dispersion (~variance)
# to avoid lowly varied genes from correlation based analysis

de18105jh <- apply(e18105jh,1,function(i) ((quantile(i,0.75)-quantile(i,0.25))/
                                             (quantile(i,0.75)+quantile(i,0.25))))
lowvarthr = 0.05

# png('GSE18105_var.png')
# par(mfrow=c(1,1))
# hist(de18105jh,breaks=40,xlim=c(0,0.6),xlab="Coefficient
#      of variation",cex.axis=1.25,cex.lab=1.25)
# dev.off()
e18105jhv <- e18105jh[de18105jh>=lowvarthr,]
print(paste(nrow(e18105jhv),' genes are selected for GSE18105'))

# for e18105, samples (columns) 1-98  are normal, 99-196 are tumor
# finally, calculate correlations based on the real and random gene sets 
library(Hmisc)
library(corpcor)
###############################################################
# Calculation of weighted degree values for the cancer networks
# self correlations are set to 0
# positive and negative correlations are considered separately
# connection requires spearman correlations p value <= pvthr
pvthr = 0.05

#TRYING SUBSETS !!!!!!
norposdegmeds2 <- c(); norposdegfqs2 <- c(); norposdegtqs2 <- c()
tumposdegmeds2 <- c(); tumposdegfqs2 <- c(); tumposdegtqs2 <- c()
nornegdegmeds2 <- c(); nornegdegfqs2 <- c(); nornegdegtqs2 <- c()
tumnegdegmeds2 <- c(); tumnegdegfqs2 <- c(); tumnegdegtqs2 <- c()

for (samplesize in 5:98){
    print(samplesize)
    norcarexp <- e18105jhv[,1:samplesize]; tumcarexp <- e18105jhv[,99:(98+samplesize)]
    
    # # take random sample sets of size 17
    # for (jj in 1:1000){
    #     print(jj)
    #     normalset <- sample(1:32,17)
    #     tumorset <- normalset + 32
    #     norcarexp <- e18105jhv[,normalset]; tumcarexp <- e18105jhv[,tumorset]
    
    # # test the original sample sets    
    # for (jj in 1:1){
    #     norcarexp <- e18105jhv[,1:32]; tumcarexp <- e18105jhv[,33:64]
    
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

# plot degree medians along the sample size range from 5 to 98
png(file='GSE18105_samplesizes_effect.png',width = 760, height = 480)
par(mfrow=c(1,2))
maxtqs <- c(max(norposdegtqs2),max(tumposdegtqs2)); maxboth <- max(maxtqs)
plot(norposdegmeds2,type='l',col='blue',lwd=4,ylim=c(0,maxboth+10),main='positive',
     xlim=c(5,98+10),ylab='degree values',xlab='sample size',cex=1.5)
lines(norposdegfqs2,lwd=0.5,lty=3,col='lightblue')
lines(norposdegtqs2,lwd=0.5,lty=3,col='lightblue')
lines(tumposdegmeds2,lwd=4,col='red')
lines(tumposdegfqs2,lwd=0.5,lty=3,col='pink')
lines(tumposdegtqs2,lwd=0.5,lty=3,col='pink')
legend('topright',legend=c('normal','tumor'),fill=c('blue','red'))

maxtqs <- c(max(nornegdegtqs2),max(tumnegdegtqs2)); maxboth <- max(maxtqs)
plot(nornegdegmeds2,type='l',col='blue',lwd=4,ylim=c(0,maxboth+10),main='negative',
     xlim=c(5,98+10),ylab='degree values',xlab='sample size',cex=1.5)
lines(nornegdegfqs2,lwd=0.5,lty=3,col='lightblue')
lines(nornegdegtqs2,lwd=0.5,lty=3,col='lightblue')
lines(tumnegdegmeds2,lwd=4,col='red')
lines(tumnegdegfqs2,lwd=0.5,lty=3,col='pink')
lines(tumnegdegtqs2,lwd=0.5,lty=3,col='pink')
legend('topright',legend=c('normal','tumor'),fill=c('blue','red'))
dev.off()

png(file='GSE18105_degs_jetset_maxvar.png',width = 760, height = 480)
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

png(file='GSE18105_degs_expthrs_4_5_6.png',width = 760, height = 480)
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

png(file='GSE18105_degs_varthrs_04_05_06.png',width = 760, height = 480)
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

png(file='GSE18105_degs_spearman_pearson.png',width = 760, height = 480)
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

png(file='GSE18105_degs_Spearmanthrs_001_005_01.png',width = 760, height = 480)
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
png(file='GSE18105_degs.png')
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
png(file='GSE18105_adecar_pvaldiff.png',width = 800, height = 800)
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
expfname = 'GSE18105_adenoma_carcinoma_pospairs_verysig.txt'
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
            if (any(genesnoname == rownames(ncorcar$P)[j]) == TRUE){
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


