
# README
# You don't need to run the code from the beginning—just skip directly to line 191:
# load("model.RData")


library(glmnet)
library(Hmisc)
library(survival)
library(survcomp)
library(survivalROC)
set.seed(2021)

ggplotROC1 <- function(dat, breaks = c(1, 2, 3), color = mycolor, smooth = TRUE) {
  library(timeROC)
  library(ggplot2)
  
  colnames(dat) <- c("OS.time", "OS", "RiskScore")
  
  # 将时间转换成年份
  if (max(dat$OS.time) > 365) {
    dat$OS.time <- dat$OS.time / 365
  } else if (max(dat$OS.time) > 24) {
    dat$OS.time <- dat$OS.time / 12
  }
  
  ROC_dat <- data.frame()
  for (br in breaks) {
    dat.test <- dat
    dat.test$OS[dat.test$OS.time >= br] <- 0
    dat.test$OS[dat.test$OS.time < br & dat.test$OS == 0] <- NA
    dat.test <- na.omit(dat.test)
    
    if (smooth) {
      roc.res <- roc(OS ~ RiskScore, data = dat.test, smooth = smooth)
    } else {
      roc.res <- roc(OS ~ RiskScore, data = dat.test)
    }
    
    test_dat <- data.frame(
      TP = roc.res$sensitivities,
      FP = 1 - roc.res$specificities,
      Times = paste0(br, " Years AUC=", round(roc.res$auc, 2))
    )
    
    ROC_dat <- rbind(ROC_dat, test_dat)
  }
  
  ggplot(ROC_dat, aes(x = FP, y = TP, fill = Times)) +
    geom_line(aes(colour = Times), lwd = 0.75) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("False positive rate") +
    ylab("True positive rate") +
    scale_colour_manual(values = color) +
    theme_bw() +
    theme(
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = NA, colour = NA)
    )
}


ggsurvplotKM<-function(dat, title = 'Groups', 
                       lables = c(), col = mycolor,
                       risk.table = TRUE,
                       tables.height = 0.25) {
  # dat：数据框，行为样本，列为时间、状态以及分组
  # the color palette to be used. Allowed values include "hue" for the default hue color scale; "grey" for grey color palettes; brewer palettes e.g. "RdBu", "Blues", ...; or custom color palette e.g. c("blue", "red"); and scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty". See details section for more information. Can be also a numeric vector of length(groups); in this case a basic color palette is created using the function palette.
  library(ggplot2)
  library(survminer)
  library(survival)
  library(ggthemes)
  library(ggplotify)
  library(ggsci)
  
  colnames(dat) <- c("OS.time", "OS", "Groups")
  # 将时间转换成年份
  if (max(dat[, 1]) > 365) {
    dat[, 1] <- dat[, 1] / 365
  } else if (max(dat[, 1]) > 24) {
    dat[, 1] <- dat[, 1] / 12
  } else {
    dat <- dat
  }
  fit <- survfit(Surv(OS.time, OS) ~ Groups,data=dat)
  surv.fit <- ggsurvplot(fit, data = dat, palette = col,
                         pval = TRUE, 
                         pval.method = T,
                         pval.method.size = 4,
                         pval.method.coord = c(0, 0.15),
                         surv.median.line='hv',
                         linetype = 1,  ###实线虚线
                         #line_size = 20,
                         pval.coord=c(0, 0.05), 
                         pval.size = 4,
                         risk.table = risk.table,
                         risk.table.y.text = FALSE,
                         legend.title = title,
                         legend.labs = lables,
                         xlab = 'Time(years)',
                         ggtheme=theme_bw(),
                         tables.height = tables.height)
  # 将图形转换为 ggplot 对象
  if (risk.table) {
    surv.fit1 <- surv.fit$plot + 
      theme(legend.position=c(1,1), 
            legend.justification=c(1,1),
            plot.margin=unit(c(0.1, 0.15, 0, 0.15), "inches"),
            legend.background = element_rect(fill = NA, colour = NA)
            # ,
            # axis.text.x=element_blank(),
            # axis.title.x=element_blank()
      )
    
    surv.fit2 <- surv.fit$table + 
      theme(plot.title=element_blank(),
            plot.margin=unit(c(0, 0.15, 0, 0.15), "inches")) +
      ylab('')
    surv.fit <- ggpubr::ggarrange(surv.fit1,
                                  surv.fit2, 
                                  ncol = 1, 
                                  nrow = 2,
                                  heights = c(1 - tables.height, 
                                              tables.height),
                                  align = "hv")
  } else {
    surv.fit <- surv.fit$plot + 
      theme(legend.position=c(1,1), 
            legend.justification=c(1,1),
            # plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
            legend.background = element_rect(fill = NA, colour = NA))
  }
  return(surv.fit)
}



# You don't need to run the code from the beginning—just skip directly to line 191:
# load("model.RData")
##################***Training***######################
load("Trainrt.RData")
TCGA <- features_train2
gene = read.csv("bootGeneselect.csv")
gene = gene$gene
ss = intersect(gene,colnames(TCGA))
TCGA <- TCGA[,ss]
#TCGA <- t(TCGA)
cli <- read.delim("tcga_cli.txt",check.names = F,row.names = 1,header = T)

row.names(cli) = substr(row.names(cli),start = 1,stop = 12)
ss = intersect(row.names(cli),row.names(TCGA))

TCGA = TCGA[ss,]
cli = cli[ss,]

colnames(cli)
cli = cli[,c("OS","OS.time")]
colnames(cli)[1] = "status"
colnames(cli)[2] = "time"
colnames(cli)
cli$time = cli$time/365
cli <- cli[,c("time","status")]

rt = merge(cli,TCGA,by = "row.names")
row.names(rt) = rt$Row.names
rt = rt[,-1]

x.train=rt

#############ridge cox###########
# 
newfit = cv.glmnet(as.matrix(x.train[,3:ncol(x.train)]),as.matrix(x.train[,1:2]), family = "cox",alpha=0,nfolds = 5)

riskScore = predict(newfit,type='link',
                    newx=as.matrix(x.train[,3:ncol(x.train)]),
                    s = newfit$lambda.min)

###############################
outCol=colnames(rt)
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
trainRiskOut=cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk))
colnames(trainRiskOut)=c("id",outCol,"riskScore","risk")
write.table(trainRiskOut,file="riskTrain.txt",sep="\t",quote=F,row.names=F)

save(newfit,file = "model.RData")

########################################################
load("model.RData")

## Validate Model Performance
dir.create("reslut")

## KM: Train set
KM <- ggsurvplotKM(trainRiskOut[, c("time", "status", "risk")],
                   lables = c('High', 'Low'),col = c("#BB0021FF","#008280FF"),
                   title = "Train")
KM

#ggsave('reslut/TCGAKM.pdf',KM,height = 4,width = 4.5)

library(pROC)
ROC <- ggplotROC1(trainRiskOut[, c("time", "status", "riskScore")],
                  color = c("#BB0021FF", "#9E5D97", "#008280FF", "purple", "orange"),
                  breaks = c(1, 2, 3, 4, 5), smooth = T)

ROC
#ggsave('reslut/TCGAROC.pdf',ROC,height = 4.2,width = 4.2)



##################**Test**###################### 
load("Testrt.RData")
EMAT <- features_vali
EMAT <- EMAT[,gene]
EMATcli <- read.delim("tcga_cli.txt",check.names = F,row.names = 1,header = T)
row.names(EMATcli) = substr(row.names(EMATcli),start = 1,stop = 12)
ss = intersect(row.names(EMATcli),row.names(EMAT))

EMAT = EMAT[ss,]
EMATcli = EMATcli[ss,]
colnames(EMATcli)

EMATcli = EMATcli[,c("OS","OS.time")]
colnames(EMATcli)[1] = "status"
colnames(EMATcli)[2] = "time"
colnames(EMATcli)
EMATcli$time = EMATcli$time/365
EMATcli <- EMATcli[,c("time","status")]
colnames(EMATcli)
EMATcli = EMATcli[,c("time","status")]
colnames(EMATcli)[1:2] = c("time","status")

rtEMAT = merge(EMATcli,EMAT,by = "row.names")
row.names(rtEMAT) = rtEMAT$Row.names
rtEMAT = rtEMAT[,-1]

x.test=rtEMAT

##################ridge
riskScore = predict(newfit,type='link',
                    newx=as.matrix(x.test[,3:ncol(x.test)]),
                    s = newfit$lambda.min)

outCol=colnames(rtEMAT)
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
TestRiskOut=cbind(id=rownames(cbind(rtEMAT[,outCol],riskScore,risk)),cbind(rtEMAT[,outCol],riskScore,risk))
colnames(TestRiskOut)=c("id",outCol,"riskScore","risk")
write.table(TestRiskOut,file="riskEMAT.txt",sep="\t",quote=F,row.names=F)

mycolor <- ggsci::pal_npg(palette = c("nrc"), alpha =1)(8)
colnames(TestRiskOut)
KM <- ggsurvplotKM(TestRiskOut[, c("time", "status", "risk")],
                   lables = c('High', 'Low'),col = c("#BB0021FF","#008280FF"),
                   title = "Test")
KM
#ggsave('reslut/TestTKM.pdf',KM,height = 4,width = 4.5)

colnames(TestRiskOut)
ROC <- ggplotROC1(TestRiskOut[, c("time", "status", "riskScore")],
                  color = c("#BB0021FF","#9E5D97","#008280FF","purple","orange"),
                  breaks = c(1,2,3,4,5),smooth = F)

ROC
#ggsave('reslut/TestROC.pdf',ROC,height = 4.2,width = 4.2)

