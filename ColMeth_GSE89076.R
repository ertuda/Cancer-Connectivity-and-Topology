# neighborhood of mutated genes in colorectal tumors
# 2: GSE89076

randomno = 1000

setwd('C:/Users/Dell/Documents/Projects/Colorectal Topology/')

x <- read.table('CosmicCompleteDifferentialMethylation.tsv',sep='\t',header=TRUE)
xc = x[x$SITE_SUBTYPE_1=='colon',]
xch = xc[xc$METHYLATION == 'H',]
xchg = unique(xch$GENE_NAME)
xcl = xc[xc$METHYLATION == 'L',]
xclg = unique(xcl$GENE_NAME)
xchgg = setdiff(xchg,xclg)
xclgg = setdiff(xclg,xchg)

# choose lowly or highly methylated genes here
xthis = xchgg 

na <- read.table('cosmic_genenames_hgnc.txt',sep='\t',header=TRUE)
nlist = intersect(na$name,xthis)
cmg = c()
for (g in nlist){
  cmg <- c(cmg,na[na$name==g,]$entrez_id)
}

# adenoma normal
print('adenoma normal')
norpairs = read.table('GSE89076_Adenoma_pos_normal_pairs.txt',header=TRUE)
norcmg = intersect(cmg,union(norpairs$Gene1,norpairs$Gene2))
print(length(norcmg))
cmgnn = norcmg
for (g in norcmg){
  cmgnn <- c(cmgnn,norpairs[norpairs$Gene1==g,]$Gene2,norpairs[norpairs$Gene2==g,]$Gene1)
}
cmgconcountan = length(cmgnn)-length(norcmg)
cmgneighcountan = length(unique(cmgnn))

# randomize to check if the connectivity is expected randomly
randconcountan = c()
randneighcountan = c()
for (i in 1:randomno){
  randnorcmg = sample(union(norpairs$Gene1,norpairs$Gene2),length(norcmg))
  cmgnn = randnorcmg
  for (g in randnorcmg){
    cmgnn <- c(cmgnn,norpairs[norpairs$Gene1==g,]$Gene2,norpairs[norpairs$Gene2==g,]$Gene1)
  }
  randconcountan = c(randconcountan,(length(cmgnn)-length(norcmg)))
  randneighcountan = c(randneighcountan,length(unique(cmgnn)))
  #print(i)
}
loweran = length(randconcountan[randconcountan<cmgconcountan])/randomno*100
loweran2 = length(randneighcountan[randneighcountan<cmgneighcountan])/randomno*100

# adenoma tumor
print('adenoma tumor')
tumpairs = read.table('GSE89076_Adenoma_pos_tumor_pairs.txt',header=TRUE)
tumcmg = intersect(cmg,union(tumpairs$Gene1,tumpairs$Gene2))
print(length(tumcmg))
cmgtn = tumcmg
for (g in tumcmg){
  cmgtn <- c(cmgtn,tumpairs[tumpairs$Gene1==g,]$Gene2,tumpairs[tumpairs$Gene2==g,]$Gene1)
}
cmgconcountat = length(cmgtn)-length(tumcmg)
cmgneighcountat <- length(unique(cmgtn))

# randomize to check if the connectivity is expected randomly
randconcountat = c()
randneighcountat = c()
for (i in 1:randomno){
  randtumcmg = sample(union(tumpairs$Gene1,tumpairs$Gene2),length(tumcmg))
  cmgtn = randtumcmg
  for (g in randtumcmg){
    cmgtn <- c(cmgtn,tumpairs[tumpairs$Gene1==g,]$Gene2,tumpairs[tumpairs$Gene2==g,]$Gene1)
  }
  randconcountat = c(randconcountat,(length(cmgtn)-length(tumcmg)))
  randneighcountat = c(randneighcountat,length(unique(cmgtn)))
  #print(i)
}
lowerat = length(randconcountat[randconcountat<cmgconcountat])/randomno*100
lowerat2 = length(randneighcountat[randneighcountat<cmgneighcountat])/randomno*100


# carcinoma normal
print('carcinoma normal')
norpairs = read.table('GSE89076_Carcinoma_pos_normal_pairs.txt',header=TRUE)
norcmg = intersect(cmg,union(norpairs$Gene1,norpairs$Gene2))
print(length(norcmg))
cmgnn = norcmg
for (g in norcmg){
  cmgnn <- c(cmgnn,norpairs[norpairs$Gene1==g,]$Gene2,norpairs[norpairs$Gene2==g,]$Gene1)
}
cmgconcountcn = length(cmgnn)-length(norcmg)
cmgneighcountcn = length(unique(cmgnn))

# randomize to check if the connectivity is expected randomly
randconcountcn = c()
randneighcountcn = c()
for (i in 1:randomno){
  randnorcmg = sample(union(norpairs$Gene1,norpairs$Gene2),length(norcmg))
  cmgnn = randnorcmg
  for (g in randnorcmg){
    cmgnn <- c(cmgnn,norpairs[norpairs$Gene1==g,]$Gene2,norpairs[norpairs$Gene2==g,]$Gene1)
  }
  randconcountcn = c(randconcountcn,(length(cmgnn)-length(norcmg)))
  randneighcountcn = c(randneighcountcn,length(unique(cmgnn)))
  #print(i)
}
lowercn = length(randconcountcn[randconcountcn<cmgconcountcn])/randomno*100
lowercn2 = length(randneighcountcn[randneighcountcn<cmgneighcountcn])/randomno*100


# carcinoma tumor
print('carcinoma tumor')
tumpairs = read.table('GSE89076_Carcinoma_pos_tumor_pairs.txt',header=TRUE)
tumcmg = intersect(cmg,union(tumpairs$Gene1,tumpairs$Gene2))
print(length(tumcmg))
cmgtn = tumcmg
for (g in tumcmg){
  cmgtn <- c(cmgtn,tumpairs[tumpairs$Gene1==g,]$Gene2,tumpairs[tumpairs$Gene2==g,]$Gene1)
}
cmgconcountct = length(cmgtn)-length(tumcmg)
cmgneighcountct <- length(unique(cmgtn))

# randomize to check if the connectivity is expected randomly
randconcountct = c()
randneighcountct = c()
for (i in 1:randomno){
  randtumcmg = sample(union(tumpairs$Gene1,tumpairs$Gene2),length(tumcmg))
  cmgtn = randtumcmg
  for (g in randtumcmg){
    cmgtn <- c(cmgtn,tumpairs[tumpairs$Gene1==g,]$Gene2,tumpairs[tumpairs$Gene2==g,]$Gene1)
  }
  randconcountct = c(randconcountct,(length(cmgtn)-length(tumcmg)))
  randneighcountct = c(randneighcountct,length(unique(cmgtn)))
  #print(i)
}
lowerct = length(randconcountct[randconcountct<cmgconcountct])/randomno*100
lowerct2 = length(randneighcountct[randneighcountct<cmgneighcountct])/randomno*100

save(randconcountan,file='Ade2randconsnormal_meth.RData')
save(randconcountat,file='Ade2randconstumor_meth.RData')
save(randconcountcn,file='Car2randconsnormal_meth.RData')
save(randconcountct,file='Car2randconstumor_meth.RData')

save(randneighcountan,file='Ade2randneighsnormal_meth.RData')
save(randneighcountat,file='Ade2randneighstumor_meth.RData')
save(randneighcountcn,file='Car2randneighsnormal_meth.RData')
save(randneighcountct,file='Car2randneighstumor_meth.RData')

# plots
#png(file='HighMeth_targets_Con_Neigh_AdeCar2_posCors.png',width=480,height=800)
tiff(file='HighMeth_targets_Con_Neigh_AdeCar2_posCors.tiff',width=3.35, height=6,units="in",res=600)
par(mfrow = c(2,2),mar=c(2.3,4,2.3,1))

# plot connection no values 
plot(sample(52:148,length(randconcountan),replace=TRUE),randconcountan,xlim=c(0,350), 
     ylim=c(-25,max(max(randconcountan),max(randconcountat))),
     cex=0.05,xaxt='n', main='Ade 2',ylab='number of connections',xlab='',cex.lab=1.25)
lines(c(52:148),rep(cmgconcountan,length(52:148)),lwd=5)
text(100,-23,paste0(as.character(loweran),'%'),cex=0.79)
points(sample(202:298,length(randconcountat),replace=TRUE),randconcountat,cex=0.05)
lines(c(202:298),rep(cmgconcountat,length(202:298)),lwd=5)
text(250,-23,paste0(as.character(lowerat),'%'),cex=0.79)
axis(1,at=c(100,250),labels=c('N','T'),tick=FALSE,cex.axis=1.15)

plot(sample(52:148,length(randconcountcn),replace=TRUE),randconcountcn,xlim=c(0,350), 
     ylim=c(-25,max(max(randconcountcn),max(randconcountct))),
     cex=0.05,xaxt='n', main='Car 2',ylab='number of connections',xlab='',cex.lab=1.25)
lines(c(52:148),rep(cmgconcountcn,length(52:148)),lwd=5)
text(100,-23,paste0(as.character(lowercn),'%'),cex=0.79)
points(sample(202:298,length(randconcountct),replace=TRUE),randconcountct,cex=0.05)
lines(c(202:298),rep(cmgconcountct,length(202:298)),lwd=5)
text(250,-23,paste0(as.character(lowerct),'%'),cex=0.79)
axis(1,at=c(100,250),labels=c('N','T'),tick=FALSE,cex.axis=1.15)


# plot neighbor no values 
plot(sample(52:148,length(randneighcountan),replace=TRUE),randneighcountan,xlim=c(0,350), 
     ylim=c(-25,max(max(randneighcountan),max(randneighcountat))),
     cex=0.05,xaxt='n', main='Ade 2',ylab='number of neighbors',xlab='',cex.lab=1.25)
lines(c(52:148),rep(cmgneighcountan,length(52:148)),lwd=5)
text(100,-23,paste0(as.character(loweran2),'%'),cex=0.79)
points(sample(202:298,length(randneighcountat),replace=TRUE),randneighcountat,cex=0.05)
lines(c(202:298),rep(cmgneighcountat,length(202:298)),lwd=5)
text(250,-23,paste0(as.character(lowerat2),'%'),cex=0.79)
axis(1,at=c(100,250),labels=c('N','T'),tick=FALSE,cex.axis=1.15)

plot(sample(52:148,length(randneighcountcn),replace=TRUE),randneighcountcn,xlim=c(0,350), 
     ylim=c(-25,max(max(randneighcountcn),max(randneighcountct))),
     cex=0.05,xaxt='n', main='Car 2',ylab='number of neighbors',xlab='',cex.lab=1.25)
lines(c(52:148),rep(cmgneighcountcn,length(52:148)),lwd=5)
text(100,-23,paste0(as.character(lowercn2),'%'),cex=0.79)
points(sample(202:298,length(randneighcountct),replace=TRUE),randneighcountct,cex=0.05)
lines(c(202:298),rep(cmgneighcountct,length(202:298)),lwd=5)
text(250,-23,paste0(as.character(lowerct2),'%'),cex=0.79)
axis(1,at=c(100,250),labels=c('N','T'),tick=FALSE,cex.axis=1.15)

dev.off()

