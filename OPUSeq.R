library(ggplot2)
library(ggrepel)
library(stringr)
library(reshape2)
library(zoo)

## Function to make a summary of all CSV count files in the set folder. 
## "m" needs to be defined before and is the number ## of different protocols per sample times two.
## This function assumes that files are named like this:
## KAPA files
## KAPA_r1_0_01_S2_noUMI_DCS_vars.csv
## KAPA_r1_0_01_S2_noUMI_SSCS_vars.csv
## KAPA_r1_0_01_S2_PP_raw_vars.csv
## OPUSeq files
## OPUS_r1_0_01_S2_DCS_vars.csv
## OPUS_r1_0_01_S2_noUMI_DCS_vars.csv
## OPUS_r1_0_01_S2_SSCS_vars.csv
## OPUS_r1_0_01_S2_noUMI_SSCS_vars.csv
## OPUS_r1_0_01_S2_PP_raw_vars.csv

summary_file <- function(fragment){
  files = list.files(pattern="*vars.csv")
  samples = sapply(strsplit(files,split="_vars"), function(a) a[1])
  protocol = sapply(strsplit(samples,split="_r"), function(a) a[1])
  sample = str_extract(samples,'(?<=_S)[:digit:]+(?=_)')
  method = str_extract(samples,'(?<=_)raw|DCS|SSCS')
  umi = files
  for (i in 1:length(files)){
    if (method[i]!="raw"){
      if (str_detect(files[i],"noUMI")){
        umi[i] = "noUMI"
      }
      else {
        umi[i] = "UMI"
      }
    }
    else {
      umi[i] = "noUMI"
    }
  }
  samples = as.data.frame(cbind(files,sample,umi,protocol,method))
  samples$sample = as.integer(samples$sample)
  samples = samples[order(samples$protocol,samples$method,samples$umi,samples$sample),]
  row.names(samples) = 1:length(files)
  ##Add a column describing the fragmentation protocol used (fragmentase or sonication)
  samples$fragment <- fragment
  ##Add expected frequencies
  samples$exp_freq <- rep(c(0,0.0001,0.0005,0.001,0.01),m)
  ## Extract observed variant frequencies
  for (i in 1:length(samples$files)){
    Sample <- read.csv(samples$files[i],header=T)
    Sample$freq <- Sample$var_count / (Sample$ref_count + Sample$var_count)
    ## get freq of spiked in variant at 534197
    if (length(Sample$freq[which((Sample$pos=="534197")&(Sample$ref=="C")&(Sample$var=="T"))])!=0){
      samples$obs_freq_197[i] <- Sample$freq[which((Sample$pos=="534197")&(Sample$ref=="C")&(Sample$var=="T"))]}
    else {samples$obs_freq_197[i] <- NA}
    ## get freq of spiked in variant at 534242
    if (length(Sample$freq[which((Sample$pos=="534242")&(Sample$ref=="A")&(Sample$var=="G"))])!=0){
      samples$obs_freq_242[i] <- Sample$freq[which((Sample$pos=="534242")&(Sample$ref=="A")&(Sample$var=="G"))]}
    else {samples$obs_freq_242[i] <- NA}
  }
  ## Calculate log10 of the obs and exp frequencies
  samples$exp_freq <- log10(samples$exp_freq)
  samples$obs_freq_197 <- log10(samples$obs_freq_197)
  samples$obs_freq_242 <- log10(samples$obs_freq_242)
  ##Calculate variant stats
  for (i in 1:dim(samples)[1]){
    Sample <- read.csv(samples$files[i],header=T)[-1]
    Sample <- Sample[which((Sample$chr=="chr11")&(Sample$pos>534038)&(Sample$pos<534466)),]
    Sample <- Sample[which(Sample$pos!=534332),]
    Sample <- Sample[which(Sample$pos!=534197),]
    Sample <- Sample[which(Sample$pos!=534242),]
    samples$total_bases[i] <- sum(Sample$ref_count) + sum(Sample$var_count)
    samples$var_bases[i] <- sum(Sample$var_count)
    samples$incidence[i] <- samples$var_bases[i] / samples$total_bases[i]
    samples$vars[i] <- sum(Sample$var_count > 0)
    samples$avg_cov[i] <- samples$total_bases[i] / 426
  }
  return(samples)
}

## Function to make barplots of VAF per position. y-axis is limited to 1.5%. Any variants above that frequency are therefore excluded.
bar_var_freq <- function(Sample,Name,Color){
  Sample <- as.data.frame(cbind(Sample$pos,Sample$freq))
  colnames(Sample) <- c("pos","freq")
  Sample$freq <- as.numeric(Sample$freq)
  Sample$pos <- as.numeric(Sample$pos)
  theplot <- ggplot(data=Sample, aes(x=pos, y=freq)) + 
    geom_col(color=Color, fill=Color) + 
    geom_col(data=Sample[which(((Sample$pos==534197)|(Sample$pos==534242))),], aes(x=pos, y=freq), fill="red", color="red",width=0.75) +
    ylim(c(0,1.5)) + ##only plotting variants below 2%!
    xlim(c(534038,534465)) +
    theme(plot.title = element_text(size=28, hjust=0.5, face="bold"),
          axis.title = element_text(size=24),
          axis.text = element_text(size=22),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_blank(),
          legend.position = "none"
    ) + 
    xlab("chr11 position") + 
    ylab("Variant frequency, %") +
    ggtitle(Name)
  return(theplot)
}

## Function to make scatterplots of expected vs. observed frequency
obs_vs_exp <- function(sub,colors,Name){
  theplot <- ggplot(sub, aes(x=exp_freq, y=value,group=variable)) + 
    geom_point(aes(color=variable),size=10,shape=9,stroke=2) + 
    scale_color_manual(values=colors)+
    theme_classic() +
    theme(plot.title = element_text(size=24, hjust=0.5, face="bold"),
          axis.title.x = element_text(size=28),
          axis.title.y = element_text(size=28),
          axis.text.x = element_text(size=32),
          axis.text.y = element_text(size=32),
          #legend.text = element_text(size=16),
          legend.position = "none",
          legend.title = element_blank()) +
    ggtitle(Name) +
    xlim(-4.5,-1.5) +
    ylim(-4.5,-1.5) +
    xlab("Expected frequency, log10") +
    ylab("Observed frequency, log10") +
    geom_abline(intercept = 0, slope = 1, colour="lightgrey", linetype="dashed", size=2)
  return(theplot)
}

##Plot var freq in "normal"/"background" sample vs that in "tumor"/"test" sample
tumor_vs_normal_freq <- function(Sample) {
  Sample <- Sample[Sample$freq_1<30,]
  ggplot(Sample, aes(x=freq_1, y=freq_2)) + geom_point(size=2) + theme_classic() + 
    geom_point(data=Sample[which(((Sample$pos==534197)&(Sample$var=="T"))|((Sample$pos==534242)&(Sample$var=="G"))),], 
               aes(x=freq_1, y=freq_2), color='red',size=3) + 
    theme(plot.title = element_text(size=24, hjust=0.5, face="bold"),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20)) +
    xlab("VF in normal, %") +
    ylab("VF in tumor, %")
}

##FDR per position plots
fdr_pos <- function(Sample){
  Sample <- Sample[which(!((Sample$pos==534332)&(Sample$var=="A"))),] ##remove germline var
  Sample <- Sample[which(Sample$freq_2>Sample$freq_1),] ##remove vars with tum freq below normal
  Sample$fdr <- p.adjust(Sample$chisq_p,method="fdr")
  Sample$log_fdr <- -log10(Sample$fdr)
  Sample$log_fdr[is.infinite(Sample$log_fdr)] <- 310
  ggplot(Sample, aes(x=pos, y=log_fdr)) + 
    geom_point(size=2) + 
    geom_hline(yintercept=1.3, linetype="dashed", color = "red") + 
    theme_classic() +
    geom_point(data=Sample[which(((Sample$pos==534197)&(Sample$var=="T"))|((Sample$pos==534242)&(Sample$var=="G"))),], 
               aes(x=pos, y=log_fdr), color='red', size=3) + 
    #geom_text_repel(data=subset(Sample,rank<11), aes(rank,log_p,label=pos),max.overlaps=100, size=6) +
    theme(plot.title = element_text(size=24, hjust=0.5, face="bold"),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text.x = element_text(size=20, colour="black"),
          axis.text.y = element_text(size=20, colour="black")) +
    xlim(c(534038,534465)) + 
    xlab("chr11 position") +
    ylab("-log10(adj. p-value)")
}

## Set directory and m
dir_kapa_plus = "..."
dir_opus_plus = "..."
dir_kapa_prep = "..."
dir_opus_prep = "..."

##Set directory containing the CSV tables 
setwd(dir_opus_plus)
##Set m (number of different computational protocols times 2. E.g. filtered reads + SSCS + DCS = 3 protocols, m is 6)
m = 10 
## Make summary table over samples. Use "frag" argument for HyperPlus files, "sonic" for HyperPrep. 
##Repeat for all datasets (OPUS and KAPA, HyperPlus and HyperPrep)
samples <- summary_file("frag")
opus_plus <- samples

## Frequency per position barplots
setwd(dir_opus_plus)
samples <- opus_plus
for (i in 1:dim(samples)[1]){
  Sample <- read.csv(samples$files[i],header=T)[-1]
  Sample$freq <- Sample$var_count * 100 / (Sample$ref_count + Sample$var_count)
  Name <- paste(samples$protocol[i],samples$sample[i],samples$method[i],samples$umi[i],sep="_")
  pdf(paste("barplots/",Name,"_var_freq.pdf",sep=""))
  print(bar_var_freq(Sample,Name,"black"))
  dev.off()
}

## Expected vs observed freq plots
setwd(dir_kapa_plus)
samples <- kapa_plus
colors <- c("#ff8055","#55aaff")
##Have to vary the Name and parameters manually here to get plots for each category
Name = "KAPA_raw"
sub <- melt(samples,id.vars=c("sample","umi","protocol","method","exp_freq"),measure.vars=c("obs_freq_197","obs_freq_242"))
sub <- sub[which(sub$exp_freq!=-Inf),]
sub <- subset(sub,protocol=="KAPA"&method=="raw"&umi=="noUMI")
pdf(paste("obs_vs_exp/",Name,"_obs_vs_exp.pdf",sep=""))
obs_vs_exp(sub,colors,Name)
dev.off()

##Plot of mapping/on target statistics

setwd("...")
files = list.files(pattern="*stats.txt")
stats = as.data.frame(1:length(files))

for (i in 1:length(files)){
  name = strsplit(files[i],split="_stats")[[1]][1]
  stats$method[i] = strsplit(name,split="_")[[1]][1]
  stats$sample[i] = str_extract(name,'(?<=_S)[:digit:]+')
  stats$name[i] = paste(stats$method[i],stats$sample[i],sep="_")
  Sample <- read.csv(files[i],sep=" ",header=F)
  stats$total[i] = Sample[1,2] / 1000000
  stats$mapped[i] = Sample[2,2] / 1000000
  stats$ontarget[i] = Sample[4,2] / 1000000
}

stats <- stats[,2:7]
stats$percent <- stats$ontarget*100 / stats$mapped
stats$sample <- as.integer(stats$sample)
stats <- stats[order(stats$method,stats$sample),]
stats_melt <- melt(stats,id.vars="name",measure.vars=c("total","mapped","ontarget"))
stats_melt$order <- rep(1:20,times=3)
stats_melt$name <- reorder(stats_melt$name,stats_melt$order)

colors=c("#55aaff","#aaaa00","#ff8055")

theplot = ggplot(data=stats_melt, aes(x=name,y=value,fill=variable)) + 
  geom_bar(stat="identity", width=0.6, position=position_dodge(width=0.8)) +
  scale_fill_manual(values=colors) +
  theme_classic() +
  xlab("Sample") +
  ylab("Million reads") +
  theme(plot.title = element_text(size=24, hjust=0.5, face="bold"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = "bottom",
        legend.title = element_blank())

ggsave(plot = theplot, width = 10, height = 6, dpi = 600, filename = "mapstats.pdf")

##Family size plot
palette3 = c('#00429d', '#1f69ad', '#3e90bd', '#5db7cd', '#7ddfdd', '#f5cb36', '#f89241', '#fc3f52', '#ce004e', '#93003a')
setwd("...")
files = list.files(pattern="*tagstats.txt")
files = files[which(!(grepl("noUMI", files, fixed = TRUE)))]
tstats1 <- read.csv(files[which((grepl("S1.", files, fixed = TRUE)))],sep="\t",header=F)[,c(1,3)]
tstats2 <- read.csv(files[which((grepl("S2.", files, fixed = TRUE)))],sep="\t",header=F)[,c(1,3)]
tstats3 <- read.csv(files[which((grepl("S3.", files, fixed = TRUE)))],sep="\t",header=F)[,c(1,3)]
tstats4 <- read.csv(files[which((grepl("S4.", files, fixed = TRUE)))],sep="\t",header=F)[,c(1,3)]
tstats5 <- read.csv(files[which((grepl("S5.", files, fixed = TRUE)))],sep="\t",header=F)[,c(1,3)]
tstats6 <- read.csv(files[which((grepl("S6.", files, fixed = TRUE)))],sep="\t",header=F)[,c(1,3)]
tstats7 <- read.csv(files[which((grepl("S7.", files, fixed = TRUE)))],sep="\t",header=F)[,c(1,3)]
tstats8 <- read.csv(files[which((grepl("S8.", files, fixed = TRUE)))],sep="\t",header=F)[,c(1,3)]
tstats9 <- read.csv(files[which((grepl("S9.", files, fixed = TRUE)))],sep="\t",header=F)[,c(1,3)]
tstats10 <- read.csv(files[which((grepl("S10.", files, fixed = TRUE)))],sep="\t",header=F)[,c(1,3)]

merged <- Reduce(function(x,y) merge(x = x, y = y, by = "V1", all=T), 
       list(tstats1,tstats2,tstats3,tstats4,tstats5,tstats6,
            tstats7,tstats8,tstats9,tstats10))
merged[is.na(merged)] = 0
colnames(merged) <- c("size","s1","s2","s3","s4","s5",
                      "s6","s7","s8","s9","s10")
merged <- merged[2:300,]

##Create a dataframe which contains the moving mean of family sizes
movmean <- data.frame(merged$size)
movmean$s1 <- c(merged$s1[1:3],rollmean(merged$s1, 4))
movmean$s2 <- c(merged$s2[1:3],rollmean(merged$s2, 4))
movmean$s3 <- c(merged$s3[1:3],rollmean(merged$s3, 4))
movmean$s4 <- c(merged$s4[1:3],rollmean(merged$s4, 4))
movmean$s5 <- c(merged$s5[1:3],rollmean(merged$s5, 4))
movmean$s6 <- c(merged$s6[1:3],rollmean(merged$s6, 4))
movmean$s7 <- c(merged$s7[1:3],rollmean(merged$s7, 4))
movmean$s8 <- c(merged$s8[1:3],rollmean(merged$s8, 4))
movmean$s9 <- c(merged$s9[1:3],rollmean(merged$s9, 4))
movmean$s10 <- c(merged$s10[1:3],rollmean(merged$s10, 4))

melted <- melt(merged,id.vars="size")
movmean <- melt(movmean,id.vars="merged.size")

theplot <- ggplot(data=melted,aes(x=size, y=value, group=variable)) + 
  geom_line(data=movmean,aes(x=merged.size,y=value,color=variable),size=1.5) + 
  scale_color_manual(values=palette3) +
  theme_classic() + 
  xlab("Family size") +
  ylab("Fraction of reads") +
  theme(axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

ggsave(theplot,filename="family_size.pdf",width=10,height=6,dpi=600)

#########################################################

## Analysis with somatic comparison

## Function to make a data frame comparing a "tumor" and a "normal" sample:
compare <- function(samples,normal,tumor){
  
  Sample = read.csv(samples$files[normal],header=T,row.names=NULL)[,-1]
  Sample <- Sample[which((Sample$chr=="chr11")&(Sample$pos>534038)&(Sample$pos<534466)),]
  Sample$qual_diff = Sample$ref_qual-Sample$var_qual
  Sample$freq = Sample$var_count*100 / (Sample$ref_count + Sample$var_count)
  normal = Sample
  
  Sample = read.csv(samples$files[tumor],header=T,row.names=NULL)[,-1]
  Sample <- Sample[which((Sample$chr=="chr11")&(Sample$pos>534038)&(Sample$pos<534466)),]
  Sample$qual_diff = Sample$ref_qual-Sample$var_qual
  Sample$freq = Sample$var_count*100 / (Sample$ref_count + Sample$var_count)
  tumor = Sample
  
  merged <- merge(normal,tumor,by=c("chr","pos","ref","var"),all=T)
  colnames(merged) <- c("chr","pos","ref","var","ref_count_1","ref_qual_1",
                        "var_count_1","var_qual_1","qual_diff_1","freq_1",
                        "ref_count_2","ref_qual_2","var_count_2","var_qual_2",
                        "qual_diff_2","freq_2")
  
  ## To fill in zeroes instead of NAs where appropriate.
  merged$var_count_1[is.na(merged$var_count_1)] = 0
  merged$var_count_2[is.na(merged$var_count_2)] = 0
  merged$freq_1[is.na(merged$freq_1)] = 0
  merged$freq_2[is.na(merged$freq_2)] = 0
  
  ## To fill in ref base coverage where it's missing in the merged DF. 
  for (i in 1:dim(merged)[1]){
    if (is.na(merged$ref_count_1[i])){
      index = 1:dim(merged)[1]
      index = index[(merged$pos == merged$pos[i])&!(is.na(merged$ref_count_1))][1]
      merged$ref_count_1[i] = merged$ref_count_1[index]
    }
    if (is.na(merged$ref_count_2[i])){
      index = 1:dim(merged)[1]
      index = index[(merged$pos == merged$pos[i])&!(is.na(merged$ref_count_2))][1]
      merged$ref_count_2[i] = merged$ref_count_2[index]
    }
    if (is.na(merged$ref_qual_1[i])){
      index = 1:dim(merged)[1]
      index = index[(merged$pos == merged$pos[i])&!(is.na(merged$ref_qual_1))][1]
      merged$ref_qual_1[i] = merged$ref_qual_1[index]
    }
    if (is.na(merged$ref_qual_2[i])){
      index = 1:dim(merged)[1]
      index = index[(merged$pos == merged$pos[i])&!(is.na(merged$ref_qual_2))][1]
      merged$ref_qual_2[i] = merged$ref_qual_2[index]
    }
  }
  
  merged <- merged[which(merged$var_count_2!=0),]
  
  ## Calculate p-value
  for (i in 1:dim(merged)[1]){
    DF <- as.data.frame(cbind(c(merged$ref_count_1[i],merged$ref_count_2[i])))
    DF$var <- c(merged$var_count_1[i],merged$var_count_2[i])
    merged$chisq_p[i] <- chisq.test(DF)$p.value
    DF <- NULL
  }
  
  ## Calculate -log10(p-val), rank by p-value, and calculate the ratio of normal 
  ## and tumor frequencies:
  
  merged$log_p <- -log10(merged$chisq_p)
  merged <- merged[order(merged$log_p,decreasing=T),]
  merged$rank <- 1:dim(merged)[1]
  merged$ratio_freq <- merged$freq_2 / merged$freq_1
  merged$LFC <- log2(merged$ratio_freq)
  
  return(merged)
}

## Isolate one set of 10 samples from a particular protocol/method
setwd(dir_opus_plus)
samples <- opus_plus
samples <- subset(samples,method=="raw"&umi=="noUMI")

norm_list <- c(1,1,1,1,6,6,6,6)
tum_list <- c(2,3,4,5,7,8,9,10)

for (i in 1:8){
  n = norm_list[i]
  t = tum_list[i]
  Name <- paste(samples$protocol[t],samples$sample[t],samples$method[t],samples$umi[t],sep="_")
  som <- compare(samples,n,t)
  write.csv(som,paste("fdr/",Name,"_fdr.csv",sep=""))
  theplot <- fdr_pos(som)
  ggsave(paste("fdr/",Name,"_fdr_pos.pdf",sep=""),plot=theplot,width=150,height=150,units="mm")
}

##Plot coverage stats

setwd(dir_opus_plus)
samples <- opus_plus
samples <- samples[which((samples$umi=="UMI")|(samples$method=="raw")),]
samples$order <- rep(c(3,1,2),each=10)
samples$method <- reorder(samples$method,samples$order)
samples$avg_cov <- log10(samples$avg_cov)
opus_prep_cov <- samples

colors=c("#aaaa00","#55aaff","#ff8055")

theplot = ggplot(data=samples, aes(x=sample,y=avg_cov,fill=method)) + 
  geom_bar(stat="identity", width=0.5, position=position_dodge(width=0.7)) +
  scale_fill_manual(values=colors) +
  theme_classic() +
  xlab("Sample") +
  ylab("log10(coverage)") +
  theme(plot.title = element_text(size=24, hjust=0.5, face="bold"),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.text = element_text(size=24),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_y_continuous(breaks=c(0:6),labels=c("1","10","100","1000","10000","100000","1000000")) +
  scale_x_continuous(breaks=c(1:10),labels=c("0","0.01","0.05","0.1","1",
                                             "0","0.01","0.05","0.1","1"))

ggsave(plot = theplot, width = 10, height = 6, dpi = 600, filename = "coverage.pdf")

## Plot variant incidence for the fragmentation method and DNA input test
summary_file <- function(){
  files = list.files(pattern="*.csv")
  samples = sapply(strsplit(files,split="_vars"), function(a) a[1])
  protocol = sapply(strsplit(samples,split="_S"), function(a) a[1])
  protocol = sapply(strsplit(protocol,split="_"), function(a) a[2])
  sample = str_extract(samples,'(?<=_S)[:digit:]+(?=_)')
  method = str_extract(samples,'(?<=_)raw|DCS|SSCS')
  dna = sapply(strsplit(samples,split="_p"), function(a) a[1])
  umi = files
  for (i in 1:length(files)){
    if (method[i]!="raw"){
      if (str_detect(files[i],"noUMI")){
        umi[i] = "noUMI"
      }
      else {
        umi[i] = "UMI"
      }
    }
    else {
      umi[i] = "noUMI"
    }
  }
  samples = as.data.frame(cbind(files,sample,dna,protocol,method,umi))
  samples$sample = as.integer(samples$sample)
  samples = samples[order(samples$sample,samples$method,samples$umi),]
  row.names(samples) = 1:dim(samples)[1]
  
  return(samples)
}
setwd(...)

## The file names were in the format: 
## na_plusold_S1_DCS_vars.csv 
## hs_plusold_S2_DCS_vars.csv
## na_plusnew_S3_DCS_vars.csv
## hs_plusnew_S4_DCS_vars.csv
## na_prep_S5_DCS_vars.csv
## hs_prep_S6_DCS_vars.csv

## Make summary table
samples <- summary_file()

##Calculate variant stats
for (i in 1:dim(samples)[1]){
  Sample <- read.csv(samples$files[i],header=T)[-1]
  Sample <- Sample[which((Sample$chr=="chr11")&(Sample$pos>534038)&(Sample$pos<534466)),]
  if (samples$dna[i] == "na"){
    Sample <- Sample[which(Sample$pos!=534332),]
  }
  if (samples$dna[i] == "hs"){
    Sample <- Sample[which(Sample$pos!=534197),]
  }
  samples$total_bases[i] <- sum(Sample$ref_count) + sum(Sample$var_count)
  samples$var_bases[i] <- sum(Sample$var_count)
  samples$incidence[i] <- samples$var_bases[i] / samples$total_bases[i]
  samples$vars[i] <- sum(Sample$var_count > 0)
  samples$avg_cov[i] <- samples$total_bases[i] / 426
}

##Plot the variant incidence at DCS and SSCS level with ggplot2.

samples <- samples[,c(1,5,9)]
samples$files <- c("1","1","2","2","3","3","4","4","5","5","6","6")
##Order the variables how I want them to appear in the plot
samples$order <- c(1,1,2,2,3,3,4,4,5,5,6,6)
samples$method_order <- c(2,1,2,1,2,1,2,1,2,1,2,1)
samples$files <- reorder(samples$files,samples$order)
samples$method <- reorder(samples$method,samples$method_order)

colors=c("#55aaff","#ff8055")

theplot <- ggplot(data=samples, aes(x=files,y=incidence,fill=method)) + 
  geom_bar(stat="identity", width=0.5, position=position_dodge(width=0.6)) +
  scale_fill_manual(values=colors) +
  theme_classic() +
  xlab("Condition") +
  ylab("Variant incidence") +
  theme(plot.title = element_text(size=24, hjust=0.5, face="bold"),
        axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28),
        axis.text.x = element_text(size=28),
        axis.text.y = element_text(size=28),
        legend.text = element_text(size=28),
        legend.position = "bottom",
        legend.title = element_blank())

ggsave(plot = theplot, width = 10, height = 6, dpi = 600, filename = "frag_test_var_incidence.pdf")


## Plotting variant incidence vs. bases trimmed from BAM. 
## Using bamUtil, the DCS BAM files were trimmed from either start, end, or both start+end 
## by 1 to 10 bases. The resulting trimmed BAM files were converted into pileups 
## and the reference and variant bases were counted in each pileup file. 
## The resulting CSV files were placed in subfolders "trim_0", "trim_1" etc to "trim_10". 

## Define the function that makes a summary of the CSV files in each directory
summary_file <- function(){
  files = list.files(pattern="*.csv")
  samples = sapply(strsplit(files,split="_vars"), function(a) a[1])
  protocol = sapply(strsplit(samples,split="_r"), function(a) a[1])
  sample = str_extract(samples,'(?<=_S)[:digit:]+(?=_)')
  method = str_extract(samples,'(?<=_)raw|DCS|SSCS')
  umi = files
  for (i in 1:length(files)){
    if (method[i]!="raw"){
      if (str_detect(files[i],"noUMI")){
        umi[i] = "noUMI"
      }
      else {
        umi[i] = "UMI"
      }
    }
    else {
      umi[i] = "noUMI"
    }
  }
  samples = as.data.frame(cbind(files,sample,umi,protocol,method))
  samples$sample = as.integer(samples$sample)
  samples = samples[order(samples$protocol,samples$method,samples$umi,samples$sample),]
  row.names(samples) = 1:length(files)
  samples$exp_freq <- rep(c(0,0.0001,0.0005,0.001,0.01),m)
  
  ## Extract observed variant frequencies
  for (i in 1:length(samples$files)){
    Sample <- read.csv(samples$files[i],header=T)
    Sample$freq <- Sample$var_count / (Sample$ref_count + Sample$var_count)
    ## get freq of spiked in variant at 534197
    if (length(Sample$freq[which((Sample$pos=="534197")&(Sample$ref=="C")&(Sample$var=="T"))])!=0){
      samples$obs_freq_197[i] <- Sample$freq[which((Sample$pos=="534197")&(Sample$ref=="C")&(Sample$var=="T"))]}
    else {samples$obs_freq_197[i] <- NA}
    ## get freq of spiked in variant at 534242
    if (length(Sample$freq[which((Sample$pos=="534242")&(Sample$ref=="A")&(Sample$var=="G"))])!=0){
      samples$obs_freq_242[i] <- Sample$freq[which((Sample$pos=="534242")&(Sample$ref=="A")&(Sample$var=="G"))]}
    else {samples$obs_freq_242[i] <- NA}
  }
  ## Calculate log10 of the obs and exp frequencies
  samples$exp_freq <- log10(samples$exp_freq)
  samples$obs_freq_197 <- log10(samples$obs_freq_197)
  samples$obs_freq_242 <- log10(samples$obs_freq_242)
  
  return(samples)
}

maindir = "...\\trim_start"
m = 4
summary = as.data.frame(0:10)

for (j in 1:11){
  k = j-1
  dir = paste(maindir,"\\trim_",k,sep="")
  setwd(dir)
  samples <- summary_file()
  samples <- samples[which(samples$umi=="UMI"),]
  for (l in 1:dim(samples)[1]){
    Sample <- read.csv(samples$files[l],header=T)[-1]
    Sample <- Sample[which((Sample$chr=="chr11")&(Sample$pos>534038)&(Sample$pos<534466)),]
    Sample <- Sample[which(Sample$pos!=534332),]
    Sample <- Sample[which(Sample$pos!=534197),]
    Sample <- Sample[which(Sample$pos!=534242),]
    samples$total_bases[l] <- sum(Sample$ref_count) + sum(Sample$var_count)
    samples$var_bases[l] <- sum(Sample$var_count)
    samples$incidence[l] <- samples$var_bases[l] / samples$total_bases[l]
    samples$vars[l] <- sum(Sample$var_count > 0)
    samples$avg_cov[l] <- samples$total_bases[l] / 426
  }
  write.csv(samples,paste(maindir,"\\trim_",k,"_var_stats.csv",sep=""),quote=F,row.names = F)
  summary$avg_var_bases[j] = sum(samples$var_bases) / 10
  summary$avg_vars[j] = sum(samples$vars) / 10
  summary$incidence[j] = sum(samples$var_bases) / sum(samples$total_bases)
  summary$avg_cov[j] = sum(samples$avg_cov) / 10
}
setwd(maindir)
write.csv(summary,"var_stats_summary.csv",quote=F,row.names = F)

##Make plots
setwd("...\\trim_start")
start <- read.csv("var_stats_summary.csv")
setwd("...\\trim_end")
end <- read.csv("var_stats_summary.csv")
setwd("...\\trim_both")
both <- read.csv("var_stats_summary.csv")

all <- cbind(start$X0.10,
             start$avg_vars,start$incidence,start$avg_cov,
             end$avg_vars,end$incidence,end$avg_cov,
             both$avg_vars,both$incidence,both$avg_cov)

colnames(all) <- c("bases_trimmed",
                   "start_vars","start_inc","start_cov",
                   "end_vars","end_inc","end_cov",
                   "both_vars","both_inc","both_cov")

all <- as.data.frame(all)
all_inc <- melt(all, id.vars = "bases_trimmed", measure.vars = c("start_inc", "end_inc", "both_inc"))
setwd("...")

colors = c("#ff8055","#17e8ad","#55aaff")
theplot <- ggplot(all_inc,aes(x=bases_trimmed, y=value, col=variable)) +
  geom_line(size=1) +
  geom_point(size=3) +
  scale_color_manual(values=colors) +
  theme_classic() +
  ylim(0,0.000003) +
  xlab("Bases trimmed") +
  ylab("Variant incidence") +
  theme(plot.title = element_text(size=20, hjust=0.5, face="bold"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.text = element_text(size=16),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_x_continuous(breaks=0:10)

ggsave(plot = theplot,"trim_bp_incidence.pdf",width=6,height=4)
