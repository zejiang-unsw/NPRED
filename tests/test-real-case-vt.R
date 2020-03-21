library(NPRED)
library(WASP)

#-------------------------------------------------------------------------------------
#load response and predictor variables
data(SPI.12); data(data.CI); data(Ind_AWAP.2.5)
#study grids and period
Grid = sample(Ind_AWAP.2.5,1)
Grid = 149 #c(45,117,142,149) 
SPI.12.ts <- window(SPI.12, start=c(1910,1),end=c(2009,12))
data.CI.ts <- window(data.CI, start=c(1910,1),end=c(2009,12))
#partition into two folds
folds <- cut(seq(1,nrow(SPI.12.ts)),breaks=2,labels=FALSE)
sub.cali <- which(folds==1, arr.ind=TRUE); sub.vali <- which(folds==2, arr.ind=TRUE)
#-------------------------------------------------------------------------------------
###calibration and selection
data <- list(x=SPI.12.ts[sub.cali,Grid],dp=data.CI.ts[sub.cali,])

#variance transformation - calibration
dwt <- modwt.vt(data, wf="d4", J=8, pad="zero", boundary="periodic")

#stepwise PIC selection
sel <- NPRED::stepwise.PIC(dwt$x, dwt$dp.n, method=F)
#-------------------------------------------------------------------------------------
###validation and prediction
data.val <- list(x=SPI.12.ts[sub.vali,Grid],dp=data.CI.ts[sub.vali,])

#variance transformation - validation
dwt.val <- modwt.vt.val(data.val, J=8, dwt)

#knn prediction
cpy <- sel$cpy; wt <- sel$wt
x=data$x; z=dwt$dp.n[,cpy]; zout=dwt.val$dp.n[,cpy]
mod <- knn(x, z, zout, k=5, pw=wt, extrap=T)

