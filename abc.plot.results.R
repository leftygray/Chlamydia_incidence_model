### R script for generating the key chlamydia plots of the H. Ali et al. paper: Figs 1, 2[a-d] & 3[a-b]

### this script requires manual editing by someone who understands R in order to plot results beyond 2012 [ignoring our current extrapolation to 2013]

# User specifications --------------------------------------------------------- 

# User specified output folder for processing
outputData <- "2018-12-03 09-52-47"  # User specified date folder

# Output folder - assuming correct working directory has been set
outputFolder <- file.path(getwd(),"output","figures")
dir.create(outputFolder) # Create folder if it doesn't already exist

# Specify type of image file for saved plots - pdf or eps
plotpdf <- TRUE

# Specify doing all plots 
allPlots <- TRUE

# Run main script -------------------------------------------------------------

# Load universal functions
source("code/load.library.R")

if (allPlots) {
  ### Read in model fit
  load(file.path("output",outputData,"posterior.dat"))
  
  ### Compute posterior weights
  is.weights <- dbeta(theta[,1],4000,420)/dbeta(theta[,1],207,22)*dbeta(theta[,5],1500,8)/dbeta(theta[,5],150,1)*dbeta(theta[,9],1050,10000)/dbeta(theta[,9],29,290)
  is.weights <- is.weights/sum(is.weights)
  
  ### Read in observational data
  source("code/abc.read.in.data.R") # there will be 16 warnings here; these can be ignored ... they are due to the missing data in the NNDSS test counts

  
  ## Fig 2: Raw Data
  setEPS()
  
  if (plotpdf) {
    pdf(file.path(outputFolder,"fig2.pdf"),width=6,height=3.267)
  } else {
    postscript(file.path(outputFolder,"fig2.eps"),width=6,height=3.267)
  }

  #
  par(mai=c(0.45,0.7,0.1,0.77),cex=0.7)
  
  plot(-100,-100,xlim=c(2001,(2001+nyears-1)),ylim=c(0,8.1*10^4),xlab="",ylab="",xaxt='n',yaxt='n')
  
  points(2001:(2001+nyears-1),notifications.m[2,]+notifications.m[3,]+notifications.m[4,]+notifications.m[5,]+notifications.m[6,],pch=24,cex=0.95,col="grey15",bg="grey15")
  lines(2001:(2001+nyears-1),notifications.m[2,]+notifications.m[3,]+notifications.m[4,]+notifications.m[5,]+notifications.m[6,],lwd=0.75,col="grey15")
  points(2001:(2001+nyears-1),notifications.f[2,]+notifications.f[3,]+notifications.f[4,]+notifications.f[5,]+notifications.f[6,],pch=23,cex=1.2,col="grey89",bg="grey89")
  lines(2001:(2001+nyears-1),notifications.f[2,]+notifications.f[3,]+notifications.f[4,]+notifications.f[5,]+notifications.f[6,],lwd=0.75,col="grey89")
  
  points(2001:(2001+nyears-1),(tested.m[2,]+tested.m[3,]+tested.m[4,])/10,pch=24,cex=0.95,col="grey45")
  lines(2001:(2001+nyears-1),(tested.m[2,]+tested.m[3,]+tested.m[4,])/10,lwd=0.75,col="grey45")
  points(2001:(2001+nyears-1),(tested.f[2,]+tested.f[3,]+tested.f[4,])/10,pch=23,cex=1.2,col="grey79")
  lines(2001:(2001+nyears-1),(tested.f[2,]+tested.f[3,]+tested.f[4,])/10,lwd=0.75,col="grey79")
  
  legend("topleft",c("Notifications","","","Tests","",""),pch=c(19,24,23,19,24,23),ncol=2,bty='n',cex=1,pt.cex=c(1,0.95,1.2,1,0.95,1.2),col="white",pt.bg="white")
  legend("topleft",c("","Males","Females      ","","Males","Females"),pch=c(19,24,23,19,24,23),ncol=2,bty='n',cex=1,pt.cex=c(1,0.95,1.2,1,0.95,1.2),col=c("transparent","grey15","grey89","transparent","grey45","grey79"),pt.bg=c("transparent","grey15","grey89","transparent","white","white"),adj=c(0,0.45))
  
  axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),labels=c("","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Number of chlamydia notifications",side=2,line=3.4,cex=0.9)
  axis(2,las=2,at=c(0,1,2,3,4,5,6,7,8,9)*10^4,labels=c("           0"," "," 20,000"," "," 40,000"," ","60,000","","80,000",""),tck=-0.0075,hadj=0.77,cex.axis=1,lwd.ticks=0.5)
  axis(2,las=2,at=c(0,1,2,3,4,5,6,7,8,9)*10^4,labels=c("            ","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  axis(4,las=2,at=c(0,1,2,3,4,5,6,7,8,9)*10^4,labels=c("            ","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  mtext("Number of chlamydia tests",side=4,line=3.9,cex=0.9)
  axis(4,las=2,at=c(0,1,2,3,4,5,6,7,8,9)*10^4,labels=c("0          "," "," 200,000"," "," 400,000"," ","600,000","","800,000",""),tck=-0.0075,hadj=0.25,cex.axis=1,lwd.ticks=0.5)
  
  dev.off()
  
  setEPS()
  if (plotpdf) {
    pdf(file.path(outputFolder,"fig3.pdf"),width=6,height=3.267*2)
  } else {
    postscript(file.path(outputFolder,"fig3.eps"),width=6,height=3.267*3)
  }
  
  layout(cbind(c(1,2,3)))
  
  ## Fig 3[a]: Notifications by sex & age-group: 15-24 & 25+
  par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
  
  plot(-100,-100,xlim=c(2001,(2001+nyears-1)+1),ylim=c(0,4.4*10^4),xlab="",ylab="",xaxt='n',yaxt='n')
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.not.f[,,1]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(seq(2001,(2001+nyears-1),by=0.01),rev(seq(2001,(2001+nyears-1),by=0.01))),
          c(predict(smooth.spline(2001:(2001+nyears-1),y.low),
                    seq(2001,(2001+nyears-1),by=0.01))$y,
            rev(predict(smooth.spline(2001:(2001+nyears-1),y.high),
                        seq(2001,(2001+nyears-1),by=0.01))$y)),angle=90,lwd=0.5,
          col="grey89",border="transparent",density=-1)
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
  lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
  lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1)+1,(2001+nyears-1)+1,by=0.01))$y[101],
                       predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),
        lwd=1,col="grey89",lty=1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.not.f[,,2]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=1,col="grey79",border="transparent",density=-1)
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
  lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
  lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),
                                          (2001+nyears-1)+1,by=0.01))$y[101],
                       predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),
        lwd=1,col="grey79",lty=1)
  
  y.matrix <- mock.not.m[,,1]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),
            rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),
          c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,
            rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),
          dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.not.m[,,2]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),
          angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),
            rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),
          c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,
            rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),
          dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")
  
  axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),labels=c("","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Number of notifications",side=2,line=3.4,cex=0.9)
  axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)*10^4,labels=c("           0",""," 10,000",""," 20,000",""," 30,000",""," 40,000",""," 50,000"),tck=-0.0075,hadj=0.77,cex.axis=1,lwd.ticks=0.5)
  axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)*10^4,labels=c("            ","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  axis(4,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)*10^4,labels=c("            ","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  points(2001:(2001+nyears-1),notifications.f[2,]+notifications.f[3,],pch=19,cex=0.8)
  points(2001:(2001+nyears-1),notifications.f[4,]+notifications.f[5,]+notifications.f[6,],pch=21,cex=0.9)
  points(2001:(2001+nyears-1),notifications.m[2,]+notifications.m[3,],pch=15,cex=0.8)
  points(2001:(2001+nyears-1),notifications.m[4,]+notifications.m[5,]+notifications.m[6,],pch=22,cex=0.9)
  
  ypos <- 1.29*10^4
  text(2013.75,ypos,"F 25+",cex=0.8)
  ypos <- 2.9*10^4
  text(2013.75,ypos,expression(plain("F 15\uad")*plain("24")),cex=0.8)
  ypos <- 2.13*10^4
  text(2013.75,ypos,"M 25+",cex=0.8)
  ypos <- 1.63*10^4
  text(2013.75,ypos,expression(plain("M 15\uad")*plain("24")),cex=0.8)
  
  legend("topleft","95% CIs",bty='n',cex=1)
  
  legend("topright",c(expression(plain("M 15\uad")*plain("24")),expression(plain("F 15\uad")*plain("24")),"M 25+","F 25+"),pch=c(15,19,22,21),ncol=2,bty='n',cex=1,pt.cex=c(0.8,0.8,0.9,0.9))
  legend("bottomleft",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("(2001+nyears-1)+1: extrapolated]")),bty='n',cex=0.8)
  
  mtext(expression(plain("Annual Notification Count from NNDSS Chlamydia Testing Australia\uad")*plain("Wide [Fitted]   ")),side=3,line=0.5,cex=0.9)
  box()
  
  ## Fig 3[b]: Test count by sex & age-group: 15-24 & 25+
  par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
  
  plot(-100,-100,xlim=c(2001,(2001+nyears-1)+1),ylim=c(0,6*10^5),xlab="",ylab="",xaxt='n',yaxt='n')
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.test.f[,,1]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
  lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
  lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],
                       predict(xy.high,seq(2013,(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey89",lty=1)
  
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.test.f[,,2]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=1,col="grey79",border="transparent",density=-1)
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
  lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
  lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],
                       predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),
        lwd=1,col="grey79",lty=1)
  
  
  y.matrix <- mock.test.m[,,1]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),
            rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),
          c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,
            rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),
          dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")
  
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.test.m[,,2]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),
            rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),
          c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,
            rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),
          dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")
  
  
  axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),labels=c("","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Number of tests",side=2,line=3.7,cex=0.9)
  axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)*10^5,labels=c("             0",""," 100,000",""," 200,000",""," 300,000",""," 400,000","","500,000","","600,000"),tck=-0.0075,hadj=0.8,cex.axis=1,lwd.ticks=0.5)
  axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)*10^5,labels=c("            ","","","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  axis(4,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)*10^5,labels=c("            ","","","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  points(2001:(2001+nyears-1),tested.f[2,],pch=19,cex=0.8)
  points(2001:(2001+nyears-1),tested.f[3,]+tested.f[4,],pch=21,cex=0.9)
  points(2001:(2001+nyears-1),tested.m[2,],pch=15,cex=0.8)
  points(2001:(2001+nyears-1),tested.m[3,]+tested.m[4,],pch=22,cex=0.9)
  
  ypos <-  3.58*10^5
  text(2013.75,ypos,expression(plain("F 15\uad")*plain("24")),cex=0.8)
  ypos <- 4.28*10^5
  text(2013.75,ypos,"F 25+",cex=0.8)
  ypos <- 1.94*10^5
  text(2013.75,ypos,"M 25+",cex=0.8)
  ypos <- 0.69*10^5
  text(2013.75,ypos,expression(plain("M 15\uad")*plain("24")),cex=0.8)
  
  legend("topleft","95% CIs",bty='n',cex=1)
  
  legend("topright",c(expression(plain("M 15\uad")*plain("24")),expression(plain("F 15\uad")*plain("24")),"M 25+","F 25+"),pch=c(15,19,22,21),ncol=2,bty='n',cex=1,pt.cex=c(0.8,0.8,0.9,0.9))
  legend("bottomright",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("(2001+nyears-1)+1: extrapolated]")),bty='n',cex=0.8)
  box()
  
  ## Fig 3[c]: Prevalence by sex & age-group: 15-24 & 25-29
  par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
  
  plot(-100,-100,xlim=c(2001,(2001+nyears-1)+1),ylim=c(0,0.075),xlab="",ylab="",xaxt='n',yaxt='n')
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.prev.f[,,1]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
  lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
  lines(c((2001+nyears-1),(2001+nyears-1)),
        c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101],
          predict(xy.high,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101]),lwd=1,col="grey89",lty=1)
  xx <- predict(xy.low,seq(2001,(2001+nyears-1),by=0.01))
  xy <- predict(xy.high,seq(2001,(2001+nyears-1),by=0.01))
  polygon(c(xx$x,rev(xy$x)),c(xx$y,rev(xy$y)),density=-1,col="grey89",border="grey89")
  lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq(2013,(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey89",lty=1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.prev.f[,,2]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
  lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
  lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey79",lty=1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.prev.m[,,1]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.prev.m[,,2]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")
  
  # Where prevalence data is used for validation - Original 2011 estimates
  prevalence.m.2011.given.had.sex <- c(4.7,6.6,3.7)*c(0.66,0.89,0.95)/100
  prevalence.f.2011.given.had.sex <- c(8.0,5.2,1.2)*c(0.56,0.90,0.97)/100
  prevalence.m.2011.given.had.sex.n <- c(298,527,432)
  prevalence.f.2011.given.had.sex.n <- c(742,1145,1140)
  
  # Updated 2014-2015 estimates - prevalence and 95% CI for each age group 
  # and gender. Provided by Jane Hocking. Age1 = 16-24, Age2 = 25-29.
  prevalence.m.2014.age1 <- c(4.47, 3.19, 6.23)/100
  prevalence.f.2014.age1 <- c(4.37, 3.51, 5.42)/100
  prevalence.m.2014.age2 <- c(1.10, 0.40, 2.86)/100
  prevalence.f.2014.age2 <- c(1.89, 1.22, 2.91)/100
  
  # Male prevalence 16-24
  p.obs <- (prevalence.m.2011.given.had.sex[1]*prevalence.m.2011.given.had.sex.n[1]+prevalence.m.2011.given.had.sex[2]*prevalence.m.2011.given.had.sex.n[2])/(prevalence.m.2011.given.had.sex.n[1]+prevalence.m.2011.given.had.sex.n[2])
  points(2011+0.12,p.obs,pch=15,cex=0.8)
  n <- (prevalence.m.2011.given.had.sex.n[1]+prevalence.m.2011.given.had.sex.n[2])
  low <- p.obs-2*sqrt(p.obs*(1-p.obs)/n)
  high <- p.obs+2*sqrt(p.obs*(1-p.obs)/n)
  lines(c(2011,2011)+0.12,c(low,high),lwd=0.75)
  lines(c(2011+0.025,2011-0.025)+0.12,c(low,low),lwd=0.75)
  lines(c(2011+0.025,2011-0.025)+0.12,c(high,high),lwd=0.75)
  # Add 2014-2015 data
  points(2014.5+0.12,prevalence.m.2014.age1[1],pch=15,cex=0.8)
  lines(c(2014.5,2014.5)+0.12,c(prevalence.m.2014.age1[2],prevalence.m.2014.age1[3]),lwd=0.75)
  lines(c(2014.5+0.025,2014.5-0.025)+0.12,c(prevalence.m.2014.age1[2],prevalence.m.2014.age1[2]),lwd=0.75)
  lines(c(2014.5+0.025,2014.5-0.025)+0.12,c(prevalence.m.2014.age1[3],prevalence.m.2014.age1[3]),lwd=0.75)
  
  # Female prevalence 16-24
  p.obs <- (prevalence.f.2011.given.had.sex[1]*prevalence.f.2011.given.had.sex.n[1]+prevalence.f.2011.given.had.sex[2]*prevalence.f.2011.given.had.sex.n[2])/(prevalence.f.2011.given.had.sex.n[1]+prevalence.f.2011.given.had.sex.n[2])
  points(2011,p.obs,pch=19,cex=0.8)
  n <- (prevalence.f.2011.given.had.sex.n[1]+prevalence.f.2011.given.had.sex.n[2])
  low <- p.obs-2*sqrt(p.obs*(1-p.obs)/n)
  high <- p.obs+2*sqrt(p.obs*(1-p.obs)/n)
  lines(c(2011,2011),c(low,high),lwd=0.75)
  lines(c(2011+0.025,2011-0.025),c(low,low),lwd=0.75)
  lines(c(2011+0.025,2011-0.025),c(high,high),lwd=0.75)
  # Add 2014-2015 data
  points(2014.5,prevalence.f.2014.age1[1],pch=19,cex=0.8)
  lines(c(2014.5,2014.5),c(prevalence.f.2014.age1[2],prevalence.f.2014.age1[3]),lwd=0.75)
  lines(c(2014.5+0.025,2014.5-0.025),c(prevalence.f.2014.age1[2],prevalence.f.2014.age1[2]),lwd=0.75)
  lines(c(2014.5+0.025,2014.5-0.025),c(prevalence.f.2014.age1[3],prevalence.f.2014.age1[3]),lwd=0.75)
  
  # Male prevalence 25-29
  p.obs <- prevalence.m.2011.given.had.sex[3]
  points(2011-0.12,p.obs,pch=22,cex=0.9)
  n <- prevalence.m.2011.given.had.sex.n[3]
  low <- p.obs-2*sqrt(p.obs*(1-p.obs)/n)
  high <- p.obs+2*sqrt(p.obs*(1-p.obs)/n)
  lines(c(2011,2011)-0.12,c(low,high),lwd=0.75)
  lines(c(2011+0.025,2011-0.025)-0.12,c(low,low),lwd=0.75)
  lines(c(2011+0.025,2011-0.025)-0.12,c(high,high),lwd=0.75)
  # Add 2014-2015 data
  points(2014.5-0.12,prevalence.m.2014.age2[1],pch=22,cex=0.8)
  lines(c(2014.5,2014.5)-0.12,c(prevalence.m.2014.age2[2],prevalence.m.2014.age2[3]),lwd=0.75)
  lines(c(2014.5+0.025,2014.5-0.025)-0.12,c(prevalence.m.2014.age2[2],prevalence.m.2014.age2[2]),lwd=0.75)
  lines(c(2014.5+0.025,2014.5-0.025)-0.12,c(prevalence.m.2014.age2[3],prevalence.m.2014.age2[3]),lwd=0.75)
  
  # Female prevalence 25-29
  p.obs <- prevalence.f.2011.given.had.sex[3]
  points(2011,p.obs,pch=21,cex=0.9)
  n <- prevalence.f.2011.given.had.sex.n[3]
  low <- p.obs-2*sqrt(p.obs*(1-p.obs)/n)
  high <- p.obs+2*sqrt(p.obs*(1-p.obs)/n)
  lines(c(2011,2011),c(low,high),lwd=0.75)
  lines(c(2011+0.025,2011-0.025),c(low,low),lwd=0.75)
  lines(c(2011+0.025,2011-0.025),c(high,high),lwd=0.75)
  # Add 2014-2015 data
  points(2014.5,prevalence.f.2014.age2[1],pch=21,cex=0.8)
  lines(c(2014.5,2014.5),c(prevalence.f.2014.age2[2],prevalence.f.2014.age2[3]),lwd=0.75)
  lines(c(2014.5+0.025,2014.5-0.025),c(prevalence.f.2014.age2[2],prevalence.f.2014.age2[2]),lwd=0.75)
  lines(c(2014.5+0.025,2014.5-0.025),c(prevalence.f.2014.age2[3],prevalence.f.2014.age2[3]),lwd=0.75)
  
  xticks <- seq(2001,(2001+nyears-1),by=2)
  finalYear <- (2001+nyears-1)
  axis(1,at=xticks,tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=xticks,labels=rep("",length(xticks)),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Prevalence (%)",side=2,line=2.5,cex=0.9)
  axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("     0","1.25","  2.5","3.75","     5","6.25","  7.5"),tck=-0.0075,hadj=0.55,cex.axis=1,lwd.ticks=0.5)
  axis
  axis(4,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  ypos <-  0.011
  text(finalYear+0.75,ypos,expression(plain("F 25\uad")*plain("29")),cex=0.8)
  ypos <- 0.025
  text(finalYear+0.75,ypos,expression(plain("F 15\uad")*plain("24")),cex=0.8)
  ypos <- 0.037
  text(finalYear+0.75,ypos,expression(plain("M 25\uad")*plain("29")),cex=0.8)
  ypos <- 0.056
  text(finalYear+0.75,ypos,expression(plain("M 15\uad")*plain("24")),cex=0.8)
  
  legend("topleft","95% CIs",bty='n',cex=1)
  
  #legend("top",c("                       ACCEPt: "),cex=1,bty='n')
  legend("topright",c(expression(plain("M 16\uad")*plain("24")),expression(plain("F 16\uad")*plain("24")),expression(plain("M 25\uad")*plain("29")),expression(plain("F 25\uad")*plain("29"))),pch=c(15,19,22,21),ncol=2,bty='n',cex=1,pt.cex=c(0.8,0.8,0.9,0.9))
  legend("bottomleft",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("(2001+nyears-1)+1: extrapolated]")),bty='n',cex=0.8)
  
  #legend("top",c("H06:                 "),cex=1,bty='n')
  #legend("top",c(expression(plain("F 18\uad")*plain("24"))),pch=c(8),ncol=2,bty='n',cex=1,pt.cex=c(0.9))
  box()
  
  dev.off()
  
  ## Fig 4[a]: Incidence numbers by sex & age-group: 15-24 & 25+
  
  # Create data frame for incidence results
  incidenceDF <- data.frame(year = numeric(nyears), 
                            male15_29_low = numeric(nyears),
                            male15_29_high = numeric(nyears),
                            female15_29_low = numeric(nyears),
                            female15_29_high = numeric(nyears))
  
  # Sort out incidence
  f.low <- f.high <- numeric(nyears)
  m.low <- m.high <- numeric(nyears)
  
  f.matrix <- mock.inc.f[,,4] 
  m.matrix <- mock.inc.m[,,4]
  
  for (jj in 1:nyears) {
    # Females
    olistf <- sort.list(f.matrix[,jj])
    f.low[jj] <- f.matrix[olistf[which.min(abs(cumsum(is.weights[olistf])-(1-0.95)/2))],jj]
    f.high[jj] <- f.matrix[olistf[which.min(abs(cumsum(is.weights[olistf])-(1-(1-0.95)/2)))],jj]
    
    # Males
    olistm <- sort.list(m.matrix[,jj])
    m.low[jj] <- m.matrix[olistm[which.min(abs(cumsum(is.weights[olistm])-(1-0.95)/2))],jj]
    m.high[jj] <- m.matrix[olistm[which.min(abs(cumsum(is.weights[olistm])-(1-(1-0.95)/2)))],jj]
  
  }

  # Store results in data frame and save
  incidenceDF$year <- 2001:(2001+nyears-1)
  incidenceDF$male15_29_low <- m.low
  incidenceDF$male15_29_high <- m.high
  incidenceDF$female15_29_low <- f.low
  incidenceDF$female15_29_high <- f.high
  
  # Write to file
  # write.csv(incidenceDF, file=file.path(outputFolder, 
  #           paste("Chlamydia_Incidence-", toString(2001+nyears-1), ".csv", sep = "")), 
  #           row.names = FALSE)

  write.csv(incidenceDF, file = file.path(outputFolder, "Chlamydia_Incidence.csv"), 
            row.names = FALSE)
  
  # Create plot
  setEPS()
  if (plotpdf) {
    pdf(file.path(outputFolder,"fig4.pdf"),width=12,height=3.267*2)
  } else {
    postscript(file.path(outputFolder,"fig4.eps"),width=6,height=3.267*2)
  }
  
  layout((c(1,2)))
  par(mai=c(0.45,0.90,0.10,0.10),cex=0.7)
  
  maxy <- 1.5
  
  plot(-100,-100,xlim=c(2001,(2001+nyears-1)+1),ylim=c(0,maxy*10^5),xlab="",ylab="",xaxt='n',yaxt='n')
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.inc.f[,,1] 
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
  lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
  lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey89",lty=1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.inc.f[,,2]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
  lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
  lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey79",lty=1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.inc.m[,,1]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.inc.m[,,2]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")
  
  axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),labels=c("","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Estimated annual number of incident",side=2,line=4.8,cex=0.9)
  mtext("chlamydia cases",side=2,line=3.55,cex=0.9)
  axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5)*10^5,labels=c("             0"," 50,000"," 100,000","150,000"," 200,000",""),tck=-0.0075,hadj=0.8,cex.axis=1,lwd.ticks=0.5)
  axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)*10^5,labels=c("            ","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  axis(4,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)*10^5,labels=c("            ","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  ypos <- 0.52*10^5
  text(2013.75,ypos,expression(plain("F 25+")),cex=0.8)
  ypos <- 0.655*10^5
  text(2013.75,ypos,expression(plain("F 15\uad")*plain("24")),cex=0.8)
  ypos <- 1.35*10^5
  text(2013.75,ypos,expression(plain("M 25+")),cex=0.8)
  ypos <- 0.905*10^5
  text(2013.75,ypos,expression(plain("M 15\uad")*plain("24")),cex=0.8)
  
  legend("topleft","95% CIs",bty='n',cex=1)
  
  legend("bottomleft",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("2013: extrapolated]")),bty='n',cex=0.8)
  
  
  ## Fig 4[b]: Incidence percentage by sex & age-group: 15-24 & 25+
  par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
  
  plot(-100,-100,xlim=c(2001,(2001+nyears-1)+1),ylim=c(0,0.105),xlab="",ylab="",xaxt='n',yaxt='n')
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.incper.f[,,1]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
  lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
  lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey89",lty=1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.incper.f[,,2]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
  lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
  lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey79",lty=1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.incper.m[,,1]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.incper.m[,,2]
  for (j in 1:nyears) {
    olist <- sort.list(y.matrix[,j])
    y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
    y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
  polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
  xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
  polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")
  
  axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),labels=c("","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Annual incidence (%)",side=2,line=2.4,cex=0.9)
  axis(2,las=2,at=seq(0,0.125,by=0.025),labels=c("     0","  2.5","     5","  7.5"," 10.0"," 12.5"),tck=-0.0075,hadj=0.6,cex.axis=1,lwd.ticks=0.5)
  axis(2,las=2,at=seq(0,0.125,by=0.025),labels=c("    ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  axis(4,las=2,at=seq(0,0.125,by=0.025),labels=c("  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  ypos <- 0.0018
  text(2013.75,ypos,expression(plain("F 25+")),cex=0.8)
  ypos <- 0.04325
  text(2013.75,ypos,expression(plain("F 15\uad")*plain("24")),cex=0.8)
  ypos <- 0.024
  text(2013.75,ypos,expression(plain("M 25+")),cex=0.8)
  ypos <- 0.083
  text(2013.75,ypos,expression(plain("M 15\uad")*plain("24")),cex=0.8)
  
  legend("topleft","95% CIs",bty='n',cex=1)
  
  legend("bottomleft",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("2013: extrapolated]")),bty='n',cex=0.8)
  
  dev.off()
  
}

## Fig 5: Time Independent Parameter Constraints
source("code/abc.read.in.hyperparameters.R")

setEPS()
if (plotpdf) {
  pdf(file.path(outputFolder,"fig5.pdf"),width=6,height=6)
} else {
  postscript(file.path(outputFolder,"fig5.eps"),width=6,height=6)
}


par(mai=c(0.45,0.6,0.10,0.15),cex=0.7)

layout(rbind(c(1,2,3),c(4,5,6),c(7,8,9)))

### Asymp.m rbeta(Nsim,400,40) # asymp.m [0.881,0.910,0.934]
posterior <- density(theta[,6],from=0,to=1,weights=is.weights,n=512,bw="SJ")

plot(-100,-100,xlim=c(0.86,0.95),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.86,0.95,length.out=100)
lines(xx,dbeta(xx,400,40)/max(dbeta(xx,400,40)),lwd=1,col="grey79")

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")

axis(1,at=c(0.88,0.91,0.94),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.88,0.91,0.94),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[asymp*plain(". [M]")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### Asymp.f rbeta(Nsim,400,75) # asymp.f [0.808,0.843,0.873]
posterior <- density(theta[,7],from=0,to=1,weights=is.weights,n=2048)

plot(-100,-100,xlim=c(0.79,0.89),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.79,0.89,length.out=100)
lines(xx,dbeta(xx,400,75)/max(dbeta(xx,400,75)),lwd=1,col="grey79")

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")

axis(1,at=c(0.80,0.83,0.86),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.80,0.83,0.86),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[asymp*plain(". [F]")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### Attend.symp # attend|symp [0.863,0.905,0.939] # 4000 420 [0.896,0.905,0.913]
posterior <- density(theta[,1],from=0,to=1,weights=is.weights,n=2048)

plot(-100,-100,xlim=c(0.85,0.95),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.85,0.95,length.out=100)
lines(xx,dbeta(xx,4000,420)/max(dbeta(xx,4000,420)),lwd=1,col="grey79")
lines(xx,dbeta(xx,207,22)/max(dbeta(xx,207,22)),lwd=1,col="grey79",lty=2)

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")
posterior <- density(theta[,1],from=0,to=1)
lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15",lty=2)

lines(xx,dbeta(xx,4000,320)/max(dbeta(xx,4000,320)),lwd=1,col="red",lty=3)
lines(xx,dbeta(xx,3200,420)/max(dbeta(xx,3200,420)),lwd=1,col="blue",lty=3)

axis(1,at=c(0.86,0.9,0.94),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.86,0.9,0.94),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[attend*plain("|symp.")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### prior.test.given.symp # rbeta(Nsim,83,4) # test|symp [0.901,0.957,0.987]

posterior <- density(theta[,2],from=0,to=1,weights=is.weights,n=2048)

plot(-100,-100,xlim=c(0.88,1),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.88,1,length.out=100)
lines(xx,dbeta(xx,83,4)/max(dbeta(xx,83,4)),lwd=1,col="grey79")

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")

lines(xx,dbeta(xx,500,10)/max(dbeta(xx,500,10)),lwd=1,col="red",lty=3)
lines(xx,dbeta(xx,1000,75)/max(dbeta(xx,1000,75)),lwd=1,col="blue",lty=3)

axis(1,at=c(0.91,0.95,0.99),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.91,0.95,0.99),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[test*plain("|symp.")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

## true.pos # rbeta(Nsim,110,7) # true.pos [0.891,0.943,0.975]

posterior <- density(theta[,3],from=0,to=1,weights=is.weights,n=2048)

plot(-100,-100,xlim=c(0.87,0.99),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.87,0.99,length.out=100)
lines(xx,dbeta(xx,110,7)/max(dbeta(xx,83,4)),lwd=1,col="grey79")

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")

lines(xx,dbeta(xx,410,10)/max(dbeta(xx,410,10)),lwd=1,col="red",lty=3)
lines(xx,dbeta(xx,810,70)/max(dbeta(xx,810,70)),lwd=1,col="blue",lty=3)

axis(1,at=c(0.88,0.93,0.98),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.88,0.93,0.98),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[pos.]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### false.pos rbeta(Nsim,2,250) # false.pos [0.001,0.007,0.022]

posterior <- density(theta[,4],from=0,to=1,weights=is.weights,n=2048,kernel="gaussian")

plot(-100,-100,xlim=c(0,0.025),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0,0.025,length.out=100)
lines(xx,dbeta(xx,2,250)/max(dbeta(xx,2,250)),lwd=1,col="grey79")

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")

axis(1,at=c(0,0.01,0.02),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0,0.01,0.02),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[false*plain(" pos.")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### rep rbeta(Nsim,150,1) # p.rep [0.976,0.995,0.999] # 1500 8 [0.990,0.995,0.998]
posterior <- density(theta[,5],from=0,to=1,weights=is.weights,n=2048)

plot(-100,-100,xlim=c(0.98,1),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.98,1,length.out=100)
lines(xx,dbeta(xx,1500,8)/max(dbeta(xx,1500,8)),lwd=1,col="grey79")
lines(xx,dbeta(xx,150,1)/max(dbeta(xx,150,1)),lwd=1,col="grey79",lty=2)

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")
posterior <- density(theta[,5],from=0,to=1)
lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15",lty=2)

lines(xx,dbeta(xx,800,1)/max(dbeta(xx,800,4)),lwd=1,col="red",lty=3)
lines(xx,dbeta(xx,3000,30)/max(dbeta(xx,3000,30)),lwd=1,col="blue",lty=3)

axis(1,at=c(0.98,0.99,1),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.98,0.99,1),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[rep.]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### prior.cured.after.year # rbeta(Nsim,22,26) # p.cured.after.year [0.321,0.458,0.599]

posterior <- density(theta[,8],from=0,to=1,weights=is.weights,n=2048)

plot(-100,-100,xlim=c(0.25,0.65),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.25,0.65,length.out=100)
lines(xx,dbeta(xx,22,26)/max(dbeta(xx,22,26)),lwd=1,col="grey79")

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")

axis(1,at=c(0.30,0.45,0.6),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.3,0.45,0.6),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[self*plain(" cure, n.")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### prior.cured.after.year # rbeta(Nsim,29,290) # p.background.antibiotic [0.062,0.090,0.125] # 1050 10000 [0.090,0.095,0.101]
posterior <- density(theta[,9],from=0,to=1,weights=is.weights,n=2048)

plot(-100,-100,xlim=c(0.05,0.14),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.05,0.14,length.out=100)
lines(xx,dbeta(xx,1050,10000)/max(dbeta(xx,1050,10000)),lwd=1,col="grey79")
lines(xx,dbeta(xx,29,290)/max(dbeta(xx,29,290)),lwd=1,col="grey79",lty=2)

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")
posterior <- density(theta[,9],from=0,to=1)
lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15",lty=2)

lines(xx,dbeta(xx,270,2000)/max(dbeta(xx,270,2000)),lwd=1,col="red",lty=3)
lines(xx,dbeta(xx,150,2000)/max(dbeta(xx,150,2000)),lwd=1,col="blue",lty=3)

axis(1,at=c(0.05,0.09,0.13),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.05,0.09,0.13),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[self*plain(" cure, a.")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

dev.off()

if ( 1 < 0) {
  
  if (1 < 0) {
    ## Fig 6: Incidence percentage by sex & age-group: 15-24 & 25+ ::: Alternative Priors
    pdf("~/Desktop/ChlamydiaPlots/fig6.pdf",width=6,height=3.267)
    par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
    
    plot(-100,-100,xlim=c(2001,(2001+nyears-1)),ylim=c(0,0.115),xlab="",ylab="",xaxt='n',yaxt='n')
    
    is.weights.alt <- dbeta(theta[,1],4000,420)/dbeta(theta[,1],207,22)#*dbeta(theta[,5],1500,8)/dbeta(theta[,5],150,1)*dbeta(theta[,9],1050,10000)/dbeta(theta[,9],29,290)
    is.weights.alt <- is.weights.alt/sum(is.weights.alt)
    #is.weights.alt <- (is.weights*0+1)/length(is.weights)
    
    y.low <- y.high <- numeric(nyears)
    y.matrix <- mock.incper.f[,,1]
    for (j in 1:nyears) {
      olist <- sort.list(y.matrix[,j])
      y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-0.95)/2))],j]
      y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-(1-0.95)/2)))],j]}
    polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
    xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
    xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
    lines(predict(xy.low,seq(2001,(2001+nyears-1),by=0.01)),lwd=1,col="grey89",lty=1)
    lines(predict(xy.high,seq(2001,(2001+nyears-1),by=0.01)),lwd=1,col="grey89",lty=1)
    lines(c((2001+nyears-1),(2001+nyears-1)),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101]),lwd=1,col="grey89",lty=1)
    xx <- predict(xy.low,seq(2001,(2001+nyears-1),by=0.01))
    xy <- predict(xy.high,seq(2001,(2001+nyears-1),by=0.01))
    #polygon(c(xx$x,rev(xy$x)),c(xx$y,rev(xy$y)),density=-1,col="grey89",border="grey89")
    lines(c((2001+nyears-1),(2001+nyears-1)),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101]),lwd=1,col="grey89",lty=1)
    
    y.low <- y.high <- numeric(nyears)
    y.matrix <- mock.incper.f[,,2]
    for (j in 1:nyears) {
      olist <- sort.list(y.matrix[,j])
      y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-0.95)/2))],j]
      y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-(1-0.95)/2)))],j]}
    polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
    xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
    xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
    lines(predict(xy.low,seq(2001,(2001+nyears-1),by=0.01)),lwd=1,col="grey79",lty=1)
    lines(predict(xy.high,seq(2001,(2001+nyears-1),by=0.01)),lwd=1,col="grey79",lty=1)
    lines(c((2001+nyears-1),(2001+nyears-1)),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101]),lwd=1,col="grey79",lty=1)
    
    y.low <- y.high <- numeric(nyears)
    y.matrix <- mock.incper.m[,,1]
    for (j in 1:nyears) {
      olist <- sort.list(y.matrix[,j])
      y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-0.95)/2))],j]
      y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-(1-0.95)/2)))],j]}
    polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
    xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
    xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
    polygon(c(seq((2001+nyears-1),(2001+nyears-1),by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1),by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y)),dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")
    
    y.low <- y.high <- numeric(nyears)
    y.matrix <- mock.incper.m[,,2]
    for (j in 1:nyears) {
      olist <- sort.list(y.matrix[,j])
      y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-0.95)/2))],j]
      y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-(1-0.95)/2)))],j]}
    polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
    xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
    xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
    polygon(c(seq((2001+nyears-1),(2001+nyears-1),by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1),by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y)),dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")
    
    axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
    axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),labels=c("","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
    mtext("Annual incidence (%)",side=2,line=2.4,cex=0.9)
    axis(2,las=2,at=seq(0,0.125,by=0.025),labels=c("     0","  2.5","     5","  7.5"," 10.0"," 12.5"),tck=-0.0075,hadj=0.55,cex.axis=1,lwd.ticks=0.5)
    axis(2,las=2,at=seq(0,0.125,by=0.025),labels=c("    ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
    axis(4,las=2,at=seq(0,0.125,by=0.025),labels=c("  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
    mtext("Year",side=1,line=1.65,cex=0.9)
    
    ypos <- 0.0018
    text(2013.75,ypos,expression(plain("F 25+")),cex=0.8)
    ypos <- 0.047
    text(2013.75,ypos,expression(plain("F 15\uad")*plain("24")),cex=0.8)
    ypos <- 0.023
    text(2013.75,ypos,expression(plain("M 25+")),cex=0.8)
    ypos <- 0.089
    text(2013.75,ypos,expression(plain("M 15\uad")*plain("24")),cex=0.8)
    
    legend("topleft","95% CIs",bty='n',cex=1)
    
    legend("bottomleft",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("2013: extrapolated]")),bty='n',cex=0.8)
    
    dev.off()
    
    
    
    ## Fig 7: Prevalence by sex & age-group: 15-24 & 25-29
    pdf("~/Desktop/fig7.pdf",width=6,height=3.267)
    par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
    
    plot(-100,-100,xlim=c(2001,(2001+nyears-1)),ylim=c(0,0.075),xlab="",ylab="",xaxt='n',yaxt='n')
    
    y.low <- y.high <- numeric(nyears)
    y.matrix <- mock.prev.f[,,1]
    for (j in 1:nyears) {
      olist <- sort.list(y.matrix[,j])
      y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-0.95)/2))],j]
      y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-(1-0.95)/2)))],j]}
    polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
    xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
    xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
    lines(predict(xy.low,seq(2001,(2001+nyears-1),by=0.01)),lwd=1,col="grey89",lty=1)
    lines(predict(xy.high,seq(2001,(2001+nyears-1),by=0.01)),lwd=1,col="grey89",lty=1)
    lines(c((2001+nyears-1),(2001+nyears-1)),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101]),lwd=1,col="grey89",lty=1)
    xx <- predict(xy.low,seq(2001,(2001+nyears-1),by=0.01))
    xy <- predict(xy.high,seq(2001,(2001+nyears-1),by=0.01))
    #polygon(c(xx$x,rev(xy$x)),c(xx$y,rev(xy$y)),density=-1,col="grey89",border="grey89")
    lines(c((2001+nyears-1),(2001+nyears-1)),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101]),lwd=1,col="grey89",lty=1)
    
    y.low <- y.high <- numeric(nyears)
    y.matrix <- mock.prev.f[,,2]
    for (j in 1:nyears) {
      olist <- sort.list(y.matrix[,j])
      y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-0.95)/2))],j]
      y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-(1-0.95)/2)))],j]}
    polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
    xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
    xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
    lines(predict(xy.low,seq(2001,(2001+nyears-1),by=0.01)),lwd=1,col="grey79",lty=1)
    lines(predict(xy.high,seq(2001,(2001+nyears-1),by=0.01)),lwd=1,col="grey79",lty=1)
    lines(c((2001+nyears-1),(2001+nyears-1)),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y[101]),lwd=1,col="grey79",lty=1)
    
    y.low <- y.high <- numeric(nyears)
    y.matrix <- mock.prev.m[,,1]
    for (j in 1:nyears) {
      olist <- sort.list(y.matrix[,j])
      y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-0.95)/2))],j]
      y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-(1-0.95)/2)))],j]}
    polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
    xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
    xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
    polygon(c(seq((2001+nyears-1),(2001+nyears-1),by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1),by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y)),dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")
    
    y.low <- y.high <- numeric(nyears)
    y.matrix <- mock.prev.m[,,2]
    for (j in 1:nyears) {
      olist <- sort.list(y.matrix[,j])
      y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-0.95)/2))],j]
      y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.alt[olist])-(1-(1-0.95)/2)))],j]}
    polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
    xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
    xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
    polygon(c(seq((2001+nyears-1),(2001+nyears-1),by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1),by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1),by=0.01))$y)),dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")
    
    prevalence.m.2011.given.had.sex <- c(4.7,6.6,3.7)*c(0.66,0.89,0.95)/100
    prevalence.f.2011.given.had.sex <- c(8.0,5.2,1.2)*c(0.56,0.90,0.97)/100
    prevalence.m.2011.given.had.sex.n <- c(298,527,432)
    prevalence.f.2011.given.had.sex.n <- c(742,1145,1140)
    
    p.obs <- (prevalence.m.2011.given.had.sex[1]*prevalence.m.2011.given.had.sex.n[1]+prevalence.m.2011.given.had.sex[2]*prevalence.m.2011.given.had.sex.n[2])/(prevalence.m.2011.given.had.sex.n[1]+prevalence.m.2011.given.had.sex.n[2])
    points(2011+0.12,p.obs,pch=15,cex=0.8)
    n <- (prevalence.m.2011.given.had.sex.n[1]+prevalence.m.2011.given.had.sex.n[2])
    low <- p.obs-2*sqrt(p.obs*(1-p.obs)/n)
    high <- p.obs+2*sqrt(p.obs*(1-p.obs)/n)
    lines(c(2011,2011)+0.12,c(low,high),lwd=0.75)
    lines(c(2011+0.025,2011-0.025)+0.12,c(low,low),lwd=0.75)
    lines(c(2011+0.025,2011-0.025)+0.12,c(high,high),lwd=0.75)
    
    p.obs <- (prevalence.f.2011.given.had.sex[1]*prevalence.f.2011.given.had.sex.n[1]+prevalence.f.2011.given.had.sex[2]*prevalence.f.2011.given.had.sex.n[2])/(prevalence.f.2011.given.had.sex.n[1]+prevalence.f.2011.given.had.sex.n[2])
    points(2011,p.obs,pch=19,cex=0.8)
    n <- (prevalence.f.2011.given.had.sex.n[1]+prevalence.f.2011.given.had.sex.n[2])
    low <- p.obs-2*sqrt(p.obs*(1-p.obs)/n)
    high <- p.obs+2*sqrt(p.obs*(1-p.obs)/n)
    lines(c(2011,2011),c(low,high),lwd=0.75)
    lines(c(2011+0.025,2011-0.025),c(low,low),lwd=0.75)
    lines(c(2011+0.025,2011-0.025),c(high,high),lwd=0.75)
    
    p.obs <- prevalence.m.2011.given.had.sex[3]
    points(2011-0.12,p.obs,pch=22,cex=0.9)
    n <- prevalence.m.2011.given.had.sex.n[3]
    low <- p.obs-2*sqrt(p.obs*(1-p.obs)/n)
    high <- p.obs+2*sqrt(p.obs*(1-p.obs)/n)
    lines(c(2011,2011)-0.12,c(low,high),lwd=0.75)
    lines(c(2011+0.025,2011-0.025)-0.12,c(low,low),lwd=0.75)
    lines(c(2011+0.025,2011-0.025)-0.12,c(high,high),lwd=0.75)
    
    p.obs <- prevalence.f.2011.given.had.sex[3]
    points(2011,p.obs,pch=21,cex=0.9)
    n <- prevalence.f.2011.given.had.sex.n[3]
    low <- p.obs-2*sqrt(p.obs*(1-p.obs)/n)
    high <- p.obs+2*sqrt(p.obs*(1-p.obs)/n)
    lines(c(2011,2011),c(low,high),lwd=0.75)
    lines(c(2011+0.025,2011-0.025),c(low,low),lwd=0.75)
    lines(c(2011+0.025,2011-0.025),c(high,high),lwd=0.75)
    
    #p.obs.med <- qbeta(0.5,5+1,160-5+1)
    #low <- qbeta(0.025,5+1,160-5+1)
    #high <- qbeta(0.975,5+1,160-5+1)
    #points(2003.5,p.obs.med,pch=8,cex=0.8)
    #lines(c(2003.5,2003.5)+0.12,c(low,high),lwd=0.75)
    #lines(c(2003.5+0.025,2003.5-0.025)+0.12,c(low,low),lwd=0.75)
    #lines(c(2003.5+0.025,2003.5-0.025)+0.12,c(high,high),lwd=0.75)
    
    axis(1,at=c(2001,2003,2005,2007,2009,2011,(2001+nyears-1)),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
    axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),labels=c("","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
    mtext("Prevalence (%)",side=2,line=2.5,cex=0.9)
    axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("     0","1.25","  2.5","3.75","     5","6.25","  7.5"),tck=-0.0075,hadj=0.55,cex.axis=1,lwd.ticks=0.5)
    axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
    axis(4,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
    mtext("Year",side=1,line=1.65,cex=0.9)
    
    
    legend("topleft","95% CIs",bty='n',cex=1)
    
    #legend("top",c("                       ACCEPt: "),cex=1,bty='n')
    legend("topright",c(expression(plain("M 16\uad")*plain("24")),expression(plain("F 16\uad")*plain("24")),expression(plain("M 25\uad")*plain("29")),expression(plain("F 25\uad")*plain("29"))),pch=c(15,22,19,21),ncol=2,bty='n',cex=1,pt.cex=c(0.8,0.9,0.8,0.9))
    legend("bottomleft",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("2013: extrapolated]")),bty='n',cex=0.8)
    
    #legend("top",c("H06:                 "),cex=1,bty='n')
    #legend("top",c(expression(plain("F 18\uad")*plain("24"))),pch=c(8),ncol=2,bty='n',cex=1,pt.cex=c(0.9))
    
    dev.off()
  }
  
}


## Fig 6: Prevalence and incidence under alternative priors
setEPS()
if (plotpdf) {
  pdf(file.path(outputFolder,"fig6.pdf"),width=12,height=3.267*3)
} else {
  postscript(file.path(outputFolder,"fig6.eps"),width=12,height=3.267*3)
}

layout(cbind(c(1,2,3),c(4,5,6)))

#### Original
is.weights <- dbeta(theta[,1],4000,420)/dbeta(theta[,1],207,22)*dbeta(theta[,5],1500,8)/dbeta(theta[,5],150,1)*dbeta(theta[,9],1050,10000)/dbeta(theta[,9],29,290)
is.weights <- is.weights/sum(is.weights)

par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)

plot(-100,-100,xlim=c(2001,(2001+nyears-1)+1),ylim=c(0,0.075),xlab="",ylab="",xaxt='n',yaxt='n')

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.f[,,1]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey89",lty=1)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.f[,,2]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey79",lty=1)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.m[,,1]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.m[,,2]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")

axis(1,at=c(2001,2003,2005,2007,2009,2011,(2001+nyears-1)),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),labels=c("","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
mtext("Prevalence (%)",side=2,line=2.5,cex=0.9)
axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("     0","1.25","  2.5","3.75","     5","6.25","  7.5"),tck=-0.0075,hadj=0.55,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
axis(4,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext("Year",side=1,line=1.65,cex=0.9)

ypos <-  0.011
text(2013.75,ypos,expression(plain("F 25\uad")*plain("29")),cex=0.8)
ypos <- 0.025
text(2013.75,ypos,expression(plain("F 15\uad")*plain("24")),cex=0.8)
ypos <- 0.037
text(2013.75,ypos,expression(plain("M 25\uad")*plain("29")),cex=0.8)
ypos <- 0.056
text(2013.75,ypos,expression(plain("M 15\uad")*plain("24")),cex=0.8)

legend("topleft","95% CIs",bty='n',cex=1)

legend("bottomleft",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("2013: extrapolated]")),bty='n',cex=0.8)
legend("top","Original priors",bty='n',cex=1.5)

#### P attend|symp.
is.weights.low <- dbeta(theta[,1],4000,320)/dbeta(theta[,1],207,22)
is.weights.low <- is.weights.low/sum(is.weights.low)
is.weights.high <- dbeta(theta[,1],3200,420)/dbeta(theta[,1],207,22)
is.weights.high <- is.weights.high/sum(is.weights.high)
is.weights <- dbeta(theta[,1],4000,420)/dbeta(theta[,1],207,22)*dbeta(theta[,5],1500,8)/dbeta(theta[,5],150,1)*dbeta(theta[,9],1050,10000)/dbeta(theta[,9],29,290)
is.weights <- is.weights/sum(is.weights)

par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)

plot(-100,-100,xlim=c(2001,(2001+nyears-1)+1),ylim=c(0,0.075),xlab="",ylab="",xaxt='n',yaxt='n')

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.f[,,1]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="grey89",density=-1)
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey89",lty=1)

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.f[,,2]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey79",lty=1)

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.m[,,1]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.m[,,2]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

axis(1,at=c(2001,2003,2005,2007,2009,2011,(2001+nyears-1)),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),labels=c("","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
mtext("Prevalence (%)",side=2,line=2.5,cex=0.9)
axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("     0","1.25","  2.5","3.75","     5","6.25","  7.5"),tck=-0.0075,hadj=0.55,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
axis(4,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext("Year",side=1,line=1.65,cex=0.9)

ypos <-  0.011
text(2013.75,ypos,expression(plain("F 25\uad")*plain("29")),cex=0.8)
ypos <- 0.025
text(2013.75,ypos,expression(plain("F 15\uad")*plain("24")),cex=0.8)
ypos <- 0.037
text(2013.75,ypos,expression(plain("M 25\uad")*plain("29")),cex=0.8)
ypos <- 0.056
text(2013.75,ypos,expression(plain("M 15\uad")*plain("24")),cex=0.8)

legend("topleft","95% CIs",bty='n',cex=1)

legend("bottomleft",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("2013: extrapolated]")),bty='n',cex=0.8)
legend("top",expression(plain("Alternative priors for ")*p[plain("attend|symp.")]),bty='n',cex=1.5)

#### P test|symp.
is.weights.low <- dbeta(theta[,2],1000,75)/dbeta(theta[,2],83,4)
is.weights.low <- is.weights.low/sum(is.weights.low)
is.weights.high <- dbeta(theta[,2],500,10)/dbeta(theta[,2],83,4)
is.weights.high <- is.weights.high/sum(is.weights.high)
is.weights <- dbeta(theta[,1],4000,420)/dbeta(theta[,1],207,22)*dbeta(theta[,5],1500,8)/dbeta(theta[,5],150,1)*dbeta(theta[,9],1050,10000)/dbeta(theta[,9],29,290)
is.weights <- is.weights/sum(is.weights)

par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)

plot(-100,-100,xlim=c(2001,(2001+nyears-1)+1),ylim=c(0,0.075),xlab="",ylab="",xaxt='n',yaxt='n')

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.f[,,1]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="grey89",density=-1)
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey89",lty=1)

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.f[,,2]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey79",lty=1)

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.m[,,1]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.m[,,2]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

axis(1,at=c(2001,2003,2005,2007,2009,2011,(2001+nyears-1)),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Prevalence (%)",side=2,line=2.5,cex=0.9)
axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("     0","1.25","  2.5","3.75","     5","6.25","  7.5"),tck=-0.0075,hadj=0.55,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
axis(4,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext("Year",side=1,line=1.65,cex=0.9)

ypos <-  0.011
text(2013.75,ypos,expression(plain("F 25\uad")*plain("29")),cex=0.8)
ypos <- 0.025
text(2013.75,ypos,expression(plain("F 15\uad")*plain("24")),cex=0.8)
ypos <- 0.037
text(2013.75,ypos,expression(plain("M 25\uad")*plain("29")),cex=0.8)
ypos <- 0.056
text(2013.75,ypos,expression(plain("M 15\uad")*plain("24")),cex=0.8)

legend("topleft","95% CIs",bty='n',cex=1)

legend("bottomleft",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("2013: extrapolated]")),bty='n',cex=0.8)
legend("top",expression(plain("Alternative priors for ")*p[plain("test|symp.")]),bty='n',cex=1.5)


#### P pos.
is.weights.low <- dbeta(theta[,3],810,80)/dbeta(theta[,3],110,7)
is.weights.low <- is.weights.low/sum(is.weights.low)
is.weights.high <- dbeta(theta[,3],410,10)/dbeta(theta[,3],110,7)
is.weights.high <- is.weights.high/sum(is.weights.high)
is.weights <- dbeta(theta[,1],4000,420)/dbeta(theta[,1],207,22)*dbeta(theta[,5],1500,8)/dbeta(theta[,5],150,1)*dbeta(theta[,9],1050,10000)/dbeta(theta[,9],29,290)
is.weights <- is.weights/sum(is.weights)

par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)

plot(-100,-100,xlim=c(2001,(2001+nyears-1)+1),ylim=c(0,0.075),xlab="",ylab="",xaxt='n',yaxt='n')

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.f[,,1]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="grey89",density=-1)
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey89",lty=1)

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.f[,,2]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey79",lty=1)

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.m[,,1]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.m[,,2]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

axis(1,at=c(2001,2003,2005,2007,2009,2011,(2001+nyears-1)),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),labels=c("","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
mtext("Prevalence (%)",side=2,line=2.5,cex=0.9)
axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("     0","1.25","  2.5","3.75","     5","6.25","  7.5"),tck=-0.0075,hadj=0.55,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
axis(4,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext("Year",side=1,line=1.65,cex=0.9)

ypos <-  0.011
text(2013.75,ypos,expression(plain("F 25\uad")*plain("29")),cex=0.8)
ypos <- 0.025
text(2013.75,ypos,expression(plain("F 15\uad")*plain("24")),cex=0.8)
ypos <- 0.037
text(2013.75,ypos,expression(plain("M 25\uad")*plain("29")),cex=0.8)
ypos <- 0.056
text(2013.75,ypos,expression(plain("M 15\uad")*plain("24")),cex=0.8)

legend("topleft","95% CIs",bty='n',cex=1)

legend("bottomleft",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("2013: extrapolated]")),bty='n',cex=0.8)
legend("top",expression(plain("Alternative priors for ")*p[plain("pos.")]),bty='n',cex=1.5)


#### P self cure, a.
is.weights.low <- dbeta(theta[,9],150,2000)/dbeta(theta[,9],29,290)
is.weights.low <- is.weights.low/sum(is.weights.low)
is.weights.high <- dbeta(theta[,9],270,2000)/dbeta(theta[,9],29,290)
is.weights.high <- is.weights.high/sum(is.weights.high)
is.weights <- dbeta(theta[,1],4000,420)/dbeta(theta[,1],207,22)*dbeta(theta[,5],1500,8)/dbeta(theta[,5],150,1)*dbeta(theta[,9],1050,10000)/dbeta(theta[,9],29,290)
is.weights <- is.weights/sum(is.weights)

par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)

plot(-100,-100,xlim=c(2001,(2001+nyears-1)+1),ylim=c(0,0.075),xlab="",ylab="",xaxt='n',yaxt='n')

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.f[,,1]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="grey89",density=-1)
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey89",lty=1)

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.f[,,2]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey79",lty=1)

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.m[,,1]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.m[,,2]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

axis(1,at=c(2001,2003,2005,2007,2009,2011,(2001+nyears-1)),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Prevalence (%)",side=2,line=2.5,cex=0.9)
axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("     0","1.25","  2.5","3.75","     5","6.25","  7.5"),tck=-0.0075,hadj=0.55,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
axis(4,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext("Year",side=1,line=1.65,cex=0.9)


ypos <-  0.011
text(2013.75,ypos,expression(plain("F 25\uad")*plain("29")),cex=0.8)
ypos <- 0.025
text(2013.75,ypos,expression(plain("F 15\uad")*plain("24")),cex=0.8)
ypos <- 0.037
text(2013.75,ypos,expression(plain("M 25\uad")*plain("29")),cex=0.8)
ypos <- 0.056
text(2013.75,ypos,expression(plain("M 15\uad")*plain("24")),cex=0.8)

legend("topleft","95% CIs",bty='n',cex=1)

legend("bottomleft",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("2013: extrapolated]")),bty='n',cex=0.8)
legend("top",expression(plain("Alternative priors for ")*p[plain("self cure, a.")]),bty='n',cex=1.5)


#### P rep.
is.weights.low <- dbeta(theta[,5],3000,30)/dbeta(theta[,5],150,1)
is.weights.low <- is.weights.low/sum(is.weights.low)
is.weights.high <- dbeta(theta[,5],1000,1)/dbeta(theta[,5],150,1)
is.weights.high <- is.weights.high/sum(is.weights.high)
is.weights <- dbeta(theta[,1],4000,420)/dbeta(theta[,1],207,22)*dbeta(theta[,5],1500,8)/dbeta(theta[,5],150,1)*dbeta(theta[,9],1050,10000)/dbeta(theta[,9],29,290)
is.weights <- is.weights/sum(is.weights)

par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)

plot(-100,-100,xlim=c(2001,(2001+nyears-1)+1),ylim=c(0,0.075),xlab="",ylab="",xaxt='n',yaxt='n')

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.f[,,1]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="grey89",density=-1)
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey89",lty=1)
lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey89",lty=1)

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.f[,,2]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
lines(predict(xy.low,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
lines(predict(xy.high,seq(2001,(2001+nyears-1)+1,by=0.01)),lwd=1,col="grey79",lty=1)
lines(c((2001+nyears-1)+1,(2001+nyears-1)+1),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101],predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y[101]),lwd=1,col="grey79",lty=1)

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.m[,,1]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=10,lwd=1.25,angle=-45,col="grey15",border="transparent",lty="13")

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

y.low <- y.high <- numeric(nyears)
y.matrix <- mock.prev.m[,,2]
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
xy.low <- smooth.spline(seq(2001,(2001+nyears-1)),y.low)
xy.high <- smooth.spline(seq(2001,(2001+nyears-1)),y.high)
polygon(c(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01),rev(seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))),c(predict(xy.low,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y,rev(predict(xy.high,seq((2001+nyears-1),(2001+nyears-1)+1,by=0.01))$y)),dens=30,lwd=0.5,angle=-45,col="grey15",border="transparent",lty="12")

for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.low[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="blue",density=0)
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.low[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-0.95)/2))],j]
  y.high[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights.high[olist])-(1-(1-0.95)/2)))],j]}
polygon(c(2001:(2001+nyears-1),(2001+nyears-1):2001),c(y.low,rev(y.high)),angle=90,lwd=0.5,col="grey89",border="red",density=0)

axis(1,at=c(2001,2003,2005,2007,2009,2011,(2001+nyears-1)),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(2001,2003,2005,2007,2009,2011,2013),labels=c("","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
mtext("Prevalence (%)",side=2,line=2.5,cex=0.9)
axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("     0","1.25","  2.5","3.75","     5","6.25","  7.5"),tck=-0.0075,hadj=0.55,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
axis(4,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext("Year",side=1,line=1.65,cex=0.9)


ypos <-  0.011
text(2013.75,ypos,expression(plain("F 25\uad")*plain("29")),cex=0.8)
ypos <- 0.025
text(2013.75,ypos,expression(plain("F 15\uad")*plain("24")),cex=0.8)
ypos <- 0.037
text(2013.75,ypos,expression(plain("M 25\uad")*plain("29")),cex=0.8)
ypos <- 0.056
text(2013.75,ypos,expression(plain("M 15\uad")*plain("24")),cex=0.8)

legend("topleft","95% CIs",bty='n',cex=1)

legend("bottomleft",expression(plain("[2001\uad")*plain("2013: model\uad")*plain("based, ")*plain("2013: extrapolated]")),bty='n',cex=0.8)
legend("top",expression(plain("Alternative priors for ")*p[plain("rep.")]),bty='n',cex=1.5)

dev.off()

y.med <- numeric(nyears)
y.matrix <- (mock.inc.m[,,1]+mock.inc.m[,,2]+mock.inc.f[,,1]+mock.inc.f[,,2])
for (j in 1:nyears) {
  olist <- sort.list(y.matrix[,j])
  y.med[j] <- y.matrix[olist[which.min(abs(cumsum(is.weights[olist])-0.5))],j]

}
