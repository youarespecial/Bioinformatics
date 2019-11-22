A_bulkData <- read.table(file = "A_bulk.txt", sep = "\t", header = T, row.names= 1)
A_bulkDiviData <- read.table(file = "A_bulk_devi.txt", sep = "\t", header = T, row.names= 1)
rawAbulk <- c(A_bulkData[1,1],A_bulkData[4,1],A_bulkData[7,1],A_bulkData[10,1],A_bulkData[13,1])




rawAbulk <- c(A_bulkDiviData[1,1],A_bulkDiviData[4,1],A_bulkDiviData[7,1],A_bulkDiviData[10,1],A_bulkDiviData[10,1])



rawA_C <- c(rawdata[1,1],rawdata[4,1],rawdata[7,1],rawdata[10,1])
rawB1_C <- c(rawdata[2,1],rawdata[5,1],rawdata[8,1],rawdata[11,1])
rawB2_C <- c(rawdata[3,1],rawdata[6,1],rawdata[9,1],rawdata[12,1])

pdf('rawdata_C_plot.pdf')
plot(time,rawA_C,type='b',ylim=c(40,750),col="DeepPink",cex=1)
par(new=TRUE)
plot(time,rawB1_C,type='b',ylim=c(40,750),col="DarkTurquoise",cex=1)
par(new=TRUE)
plot(time,rawB2_C,type='b',ylim=c(40,750),col="RosyBrown",cex=1)
legend(1,700,c('lung','liveHigh','liveLow'),col=c('DeepPink',"DarkTurquoise","RosyBrown"),text.col=c("DeepPink","DarkTurquoise","RosyBrown") )
dev.off()


rawA_adar1 <- c(rawdata[1,2],rawdata[4,2],rawdata[7,2],rawdata[10,2])
rawB1_adar1 <- c(rawdata[2,2],rawdata[5,2],rawdata[8,2],rawdata[11,2])
rawB2_adar1 <- c(rawdata[3,2],rawdata[6,2],rawdata[9,2],rawdata[12,2])

pdf('rawdata_adar1_plot.pdf')
plot(time,rawA_adar1,type='b',ylim=c(70,1450),col="DeepPink",cex=1)
par(new=TRUE)
plot(time,rawB1_adar1,type='b',ylim=c(70,1450),col="DarkTurquoise",cex=1)
par(new=TRUE)
plot(time,rawB2_adar1,type='b',ylim=c(70,1450),col="RosyBrown",cex=1)
legend(3.5,4,c('lung','liveHigh','liveLow'),col=c('DeepPink',"DarkTurquoise","RosyBrown"),text.col=c("DeepPink","DarkTurquoise","RosyBrown") )
dev.off()

rawA_adarb2 <- c(rawdata[1,3],rawdata[4,3],rawdata[7,3],rawdata[10,3])
rawB1_adarb2 <- c(rawdata[2,3],rawdata[5,3],rawdata[8,3],rawdata[11,3])
rawB2_adarb2 <- c(rawdata[3,3],rawdata[6,3],rawdata[9,3],rawdata[12,3])

pdf('rawdata_adarb2_plot.pdf')
plot(time,rawA_adarb2,type='b',ylim=c(30,550),col="DeepPink",cex=1)
par(new=TRUE)
plot(time,rawB1_adarb2,type='b',ylim=c(30,550),col="DarkTurquoise",cex=1)
par(new=TRUE)
plot(time,rawB2_adarb2,type='b',ylim=c(30,550),col="RosyBrown",cex=1)
legend(3.5,4,c('lung','liveHigh','liveLow'),col=c('DeepPink',"DarkTurquoise","RosyBrown"),text.col=c("DeepPink","DarkTurquoise","RosyBrown") )
dev.off()

rawA_adar12 <- c(rawdata[1,4],rawdata[4,4],rawdata[7,4],rawdata[10,4])
rawB1_adar12 <- c(rawdata[2,4],rawdata[5,4],rawdata[8,4],rawdata[11,4])
rawB2_adar12 <- c(rawdata[3,4],rawdata[6,4],rawdata[9,4],rawdata[12,4])

pdf('rawdata_adar12_plot.pdf')
plot(time,rawA_adar12,type='b',ylim=c(0,480),col="DeepPink",cex=1)
par(new=TRUE)
plot(time,rawB1_adar12,type='b',ylim=c(0,480),col="DarkTurquoise",cex=1)
par(new=TRUE)
plot(time,rawB2_adar12,type='b',ylim=c(0,480),col="RosyBrown",cex=1)
legend(1,480,c('lung','liveHigh','liveLow'),col=c('DeepPink',"DarkTurquoise","RosyBrown"),text.col=c("DeepPink","DarkTurquoise","RosyBrown") )
dev.off()
#####################################

pdf('A_plot.pdf')
plot(time,rawA_C,type='b',ylim=c(0,250),col="DeepPink",cex=1)
par(new=TRUE)
plot(time,rawA_adar1,type='b',ylim=c(0,250),col="RosyBrown",cex=1)
par(new=TRUE)
plot(time,rawA_adarb2,type='b',ylim=c(0,250),col="DarkTurquoise",cex=1)
par(new=TRUE)
plot(time,rawA_adar12,type='b',ylim=c(0,250),col="dimgrey",cex=1)
legend(1,1400,c('C','adar1','adarb2','adar12'),col=c('DeepPink',"RosyBrown","DarkTurquoise","dimgrey"),text.col=c("DeepPink","RosyBrown","DarkTurquoise","dimgrey") )
dev.off()



pdf('B1_plot.pdf')
plot(time,rawB1_C,type='b',ylim=c(0,1400),col="DeepPink",cex=1)
par(new=TRUE)
plot(time,rawB1_adar1,type='b',ylim=c(0,1400),col="RosyBrown",cex=1)
par(new=TRUE)
plot(time,rawB1_adarb2,type='b',ylim=c(0,1400),col="DarkTurquoise",cex=1)
par(new=TRUE)
plot(time,rawB1_adar12,type='b',ylim=c(0,1400),col="RosyBrown",cex=1)
legend(1,1400,c('C','adar1','adarb2','adar12'),col=c('DeepPink',"RosyBrown","DarkTurquoise","dimgrey"),text.col=c("DeepPink","RosyBrown","DarkTurquoise","dimgrey") )
dev.off()



pdf('B2_plot.pdf')
plot(time,rawB2_C,type='b',ylim=c(0,1400),col="DeepPink",cex=1)
par(new=TRUE)
plot(time,rawB2_adar1,type='b',ylim=c(0,1400),col="RosyBrown",cex=1)
par(new=TRUE)
plot(time,rawB2_adarb2,type='b',ylim=c(0,1400),col="DarkTurquoise",cex=1)
par(new=TRUE)
plot(time,rawB2_adar12,type='b',ylim=c(0,1400),col="dimgrey",cex=1)
legend(1,1400,c('C','adar1','adarb2','adar12'),col=c('DeepPink',"RosyBrown","DarkTurquoise","dimgrey"),text.col=c("DeepPink","RosyBrown","DarkTurquoise","dimgrey") )
dev.off()




#################ABC############
pdf('AB_plot.pdf')
plot(time,rawA_C,type='b',ylim=c(0,1400),col="DeepPink",cex=1,lty=1)
par(new=TRUE)
plot(time,rawA_adar1,type='b',ylim=c(0,1400),col="RosyBrown",cex=1,lty=1)
par(new=TRUE)
plot(time,rawA_adarb2,type='b',ylim=c(0,1400),col="DarkTurquoise",cex=1,lty=1)
par(new=TRUE)
plot(time,rawA_adar12,type='b',ylim=c(0,1400),col="dimgrey",cex=1,lty=1)

par(new=TRUE)

plot(time,rawB1_C,type='b',ylim=c(0,1400),col="DeepPink",cex=1,lty=2)
par(new=TRUE)
plot(time,rawB1_adar1,type='p',ylim=c(0,1400),col="RosyBrown",cex=1,lty=2)
par(new=TRUE)
plot(time,rawB1_adarb2,type='p',ylim=c(0,1400),col="DarkTurquoise",cex=1,lty=2)
par(new=TRUE)
plot(time,rawB1_adar12,type='p',ylim=c(0,1400),col="dimgrey",cex=1,lty=2)
par(new=TRUE)


plot(time,rawB2_C,type='c',ylim=c(0,1400),col="DeepPink",cex=1,lty=3)
par(new=TRUE)
plot(time,rawB2_adar1,type='c',ylim=c(0,1400),col="RosyBrown",cex=1,lty=3)
par(new=TRUE)
plot(time,rawB2_adarb2,type='c',ylim=c(0,1400),col="DarkTurquoise",cex=1,lty=3)
par(new=TRUE)
plot(time,rawB2_adar12,type='c',ylim=c(0,1400),col="dimgrey",cex=1,lty=3)


legend(1,1400,c('A_C','A_adar1','A_adarb2','A_adar12','B1_C','B1_adar1','B1_adarb2','B1_adar12','B2_C','B2_adar1','B2_adarb2','B2_adar12'),col=c('DeepPink',"RosyBrown","DarkTurquoise","dimgrey",'DeepPink',"RosyBrown","DarkTurquoise","dimgrey",'DeepPink',"RosyBrown","DarkTurquoise","dimgrey"),text.col=c("DeepPink","RosyBrown","DarkTurquoise","dimgrey"),lty=c(1,1,1,1,2,2,2,2,3,3,3,3))
dev.off()


#####################divide
dataA_adar1 <- c(data[1,1],data[4,1],data[7,1],data[10,1])
dataB1_adar1 <- c(data[2,1],data[5,1],data[8,1],data[11,1])
dataB2_adar1 <- c(data[3,1],data[6,1],data[9,1],data[12,1])
dataA_adarb2 <- c(data[1,2],data[4,2],data[7,2],data[10,2])
dataB1_adarb2 <- c(data[2,2],data[5,2],data[8,2],data[11,2])
dataB2_adarb2 <- c(data[3,2],data[6,2],data[9,2],data[12,2])
dataA_adar12 <- c(data[1,3],data[4,3],data[7,3],data[10,3])
dataB1_adar12 <- c(data[2,3],data[5,3],data[8,3],data[11,3])
dataB2_adar12 <- c(data[3,3],data[6,3],data[9,3],data[12,3])


pdf('AB_divided_plot.pdf')
plot(time,dataA_adar1,type='b',ylim=c(0,4),col="RosyBrown",cex=1,lty=1)
par(new=TRUE)
plot(time,dataA_adarb2,type='b',ylim=c(0,4),col="DarkTurquoise",cex=1,lty=1)
par(new=TRUE)
plot(time,dataA_adar12,type='b',ylim=c(0,4),col="dimgrey",cex=1,lty=1)

par(new=TRUE)
plot(time,dataB1_adar1,type='b',ylim=c(0,4),col="RosyBrown",cex=1,lty=2)
par(new=TRUE)
plot(time,dataB1_adarb2,type='b',ylim=c(0,4),col="DarkTurquoise",cex=1,lty=2)
par(new=TRUE)
plot(time,dataB1_adar12,type='b',ylim=c(0,4),col="dimgrey",cex=1,lty=2)
par(new=TRUE)

plot(time,dataB2_adar1,type='b',ylim=c(0,4),col="RosyBrown",cex=1,lty=3)
par(new=TRUE)
plot(time,dataB2_adarb2,type='b',ylim=c(0,4),col="DarkTurquoise",cex=1,lty=3)
par(new=TRUE)
plot(time,dataB2_adar12,type='b',ylim=c(0,4),col="dimgrey",cex=1,lty=3)


legend(1,4,c('A_adar1','A_adarb2','A_adar12','B1_adar1','B1_adarb2','B1_adar12','B2_adar1','B2_adarb2','B2_adar12'),col=c("RosyBrown","DarkTurquoise","dimgrey","RosyBrown","DarkTurquoise","dimgrey","RosyBrown","DarkTurquoise","dimgrey"),text.col=c("RosyBrown","DarkTurquoise","dimgrey"),lty=c(1,1,1,2,2,2,3,3,3))
dev.off()

####################
pdf('A_divided_plot.pdf')
plot(time,dataA_adar1,type='b',ylim=c(0,4),col="RosyBrown",cex=1)
par(new=TRUE)
plot(time,dataA_adarb2,type='b',ylim=c(0,4),col="DarkTurquoise",cex=1)
par(new=TRUE)
plot(time,dataA_adar12,type='b',ylim=c(0,4),col="dimgrey",cex=1)
legend(1,4,c('adar1','adarb2','adar12'),col=c("RosyBrown","DarkTurquoise","dimgrey"),text.col=c("RosyBrown","DarkTurquoise","dimgrey") )
dev.off()

pdf('B1_divided_plot.pdf')
plot(time,dataB1_adar1,type='b',ylim=c(0,4),col="RosyBrown",cex=1)
par(new=TRUE)
plot(time,dataB1_adarb2,type='b',ylim=c(0,4),col="DarkTurquoise",cex=1)
par(new=TRUE)
plot(time,dataB1_adar12,type='b',ylim=c(0,4),col="RosyBrown",cex=1)
legend(1,4,c('adar1','adarb2','adar12'),col=c("RosyBrown","DarkTurquoise","dimgrey"),text.col=c("RosyBrown","DarkTurquoise","dimgrey") )
dev.off()

pdf('B2_divided_plot.pdf')
plot(time,dataB2_adar1,type='b',ylim=c(0,4),col="RosyBrown",cex=1)
par(new=TRUE)
plot(time,dataB2_adarb2,type='b',ylim=c(0,4),col="DarkTurquoise",cex=1)
par(new=TRUE)
plot(time,dataB2_adar12,type='b',ylim=c(0,4),col="dimgrey",cex=1)
legend(1,4,c('adar1','adarb2','adar12'),col=c("RosyBrown","DarkTurquoise","dimgrey"),text.col=c("RosyBrown","DarkTurquoise","dimgrey") )
dev.off()
###################################      tian at last    ########################
pdf('A_plot.pdf')
plot(time,rawA_C,type='b',ylim=c(0,250),col="red",cex=1,pch=16,lwd=2)
par(new=TRUE)
plot(time,rawA_adar1,type='b',ylim=c(0,250),col="darkblue",cex=1,pch=16,lwd=2)
par(new=TRUE)
plot(time,rawA_adarb2,type='b',ylim=c(0,250),col="DarkTurquoise",cex=1,pch=16,lwd=2)
par(new=TRUE)
plot(time,rawA_adar12,type='b',ylim=c(0,250),col="yellowgreen",cex=1,pch=16,lwd=2)
legend(1,250,c('C','adar1','adarb2','adar12'),col=c('red',"darkblue","DarkTurquoise","yellowgreen"),text.col=c("red","darkblue","DarkTurquoise","yellowgreen") )
dev.off()

pdf('A_divided_plot.pdf')
plot(time,dataA_adar1,type='b',ylim=c(0,4),col="darkblue",cex=1,pch=16,lwd=2)
par(new=TRUE)
plot(time,dataA_adarb2,type='b',ylim=c(0,4),col="DarkTurquoise",cex=1,pch=16,lwd=2)
par(new=TRUE)
plot(time,dataA_adar12,type='b',ylim=c(0,4),col="red",cex=1,pch=16,lwd=2)
legend(1,4,c('adar1','adarb2','adar12'),col=c("darkblue","DarkTurquoise","red"),text.col=c("darkblue","DarkTurquoise","red"))
dev.off()
##########################相同时间，三种case相对control的t.test（）##################



> t.test(A15['C_bulk'],A15['adar1_bulk'])

	Welch Two Sample t-test

data:  A15["C_bulk"] and A15["adar1_bulk"]
t = -3.6122, df = 7.1984, p-value = 0.008199
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -237.31174  -50.16926
sample estimates:
mean of x mean of y
  99.0874  242.8279

> t.test(A15['C_bulk'],A15['adarb2_bulk'])

	Welch Two Sample t-test

data:  A15["C_bulk"] and A15["adarb2_bulk"]
t = -2.3053, df = 6.4309, p-value = 0.05779
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -215.195249    4.678249
sample estimates:
mean of x mean of y
  99.0874  204.3459

> t.test(A15['C_bulk'],A15['adar12_bulk'])

	Welch Two Sample t-test

data:  A15["C_bulk"] and A15["adar12_bulk"]
t = 0.53927, df = 6.1435, p-value = 0.6087
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -49.4145  77.5557
sample estimates:
mean of x mean of y
  99.0874   85.0168



t.test(A6['C_bulk'],A6['adar1_bulk'], var.equal=TRUE)
t.test(A6['C_bulk'],A6['adarb2_bulk'], var.equal=TRUE)
t.test(A6['C_bulk'],A6['adar12_bulk'], var.equal=TRUE)

t.test(A9['C_bulk'],A9['adar1_bulk'], var.equal=TRUE)
t.test(A9['C_bulk'],A9['adarb2_bulk'], var.equal=TRUE)
t.test(A6['C_bulk'],A9['adar12_bulk'], var.equal=TRUE)

t.test(A13['C_bulk'],A13['adar1_bulk'], var.equal=TRUE)
t.test(A13['C_bulk'],A13['adarb2_bulk'], var.equal=TRUE)
t.test(A13['C_bulk'],A13['adar12_bulk'], var.equal=TRUE)

t.test(A15['C_bulk'],A15['adar1_bulk'], var.equal=TRUE)
t.test(A15['C_bulk'],A15['adarb2_bulk'], var.equal=TRUE)
t.test(A15['C_bulk'],A15['adar12_bulk'], var.equal=TRUE)



> t.test(A6['C_bulk'],A6['adar1_bulk'])

	Welch Two Sample t-test

data:  A6["C_bulk"] and A6["adar1_bulk"]
t = -2.5561, df = 5.8075, p-value = 0.0444
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -105.360578   -1.874222
sample estimates:
mean of x mean of y
  33.4784   87.0958

> t.test(A6['C_bulk'],A6['adarb2_bulk'])

	Welch Two Sample t-test

data:  A6["C_bulk"] and A6["adarb2_bulk"]
t = -0.69303, df = 4.8266, p-value = 0.5202
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -98.58037  57.06957
sample estimates:
mean of x mean of y
  33.4784   54.2338

> t.test(A6['C_bulk'],A6['adar12_bulk'])

	Welch Two Sample t-test

data:  A6["C_bulk"] and A6["adar12_bulk"]
t = 3.6349, df = 4, p-value = 0.02206
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  7.906926 59.049874
sample estimates:
mean of x mean of y
  33.4784    0.0000

>
> t.test(A9['C_bulk'],A9['adar1_bulk'])

	Welch Two Sample t-test

data:  A9["C_bulk"] and A9["adar1_bulk"]
t = -4.2094, df = 5.2098, p-value = 0.007688
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -178.41229  -44.13751
sample estimates:
mean of x mean of y
  56.2635  167.5384

> t.test(A9['C_bulk'],A9['adarb2_bulk'])

	Welch Two Sample t-test

data:  A9["C_bulk"] and A9["adarb2_bulk"]
t = -2.6539, df = 6.269, p-value = 0.03629
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -100.864479   -4.614321
sample estimates:
mean of x mean of y
  56.2635  109.0029

> t.test(A6['C_bulk'],A9['adar12_bulk'])

	Welch Two Sample t-test

data:  A6["C_bulk"] and A9["adar12_bulk"]
t = -1.7704, df = 5.1066, p-value = 0.1356
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -42.550736   7.714136
sample estimates:
mean of x mean of y
  33.4784   50.8967

>
>
> t.test(A13['C_bulk'],A13['adar1_bulk'])

	Welch Two Sample t-test

data:  A13["C_bulk"] and A13["adar1_bulk"]
t = -3.3904, df = 5.8881, p-value = 0.0151
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -160.1627  -25.5261
sample estimates:
mean of x mean of y
  83.9521  176.7965

> t.test(A13['C_bulk'],A13['adarb2_bulk'])

	Welch Two Sample t-test

data:  A13["C_bulk"] and A13["adarb2_bulk"]
t = -2.3439, df = 6.1138, p-value = 0.05675
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -124.332532    2.393332
sample estimates:
mean of x mean of y
  83.9521  144.9217

> t.test(A13['C_bulk'],A13['adar12_bulk'])

	Welch Two Sample t-test

data:  A13["C_bulk"] and A13["adar12_bulk"]
t = 1.9831, df = 7.1266, p-value = 0.08706
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -5.568756 64.769356
sample estimates:
mean of x mean of y
  83.9521   54.3518

# 绘制箱线图
p <- ggboxplot(Adevi6, x="supp", y="len", color = "supp",
  palette = "jco", add = "jitter") 




#################### 9.02 updata all data have been detected ###########
bulkData <- read.table(file = "tumor_size.txt", sep = "\t", header = T, row.names= 1)
dataA_c <- c(bulkData[1,1],bulkData[4,1],bulkData[7,1],bulkData[10,1],bulkData[13,1])
dataA_adar1 <- c(bulkData[1,2],bulkData[4,2],bulkData[7,2],bulkData[10,2],bulkData[13,2])
dataA_adarb2 <- c(bulkData[1,3],bulkData[4,3],bulkData[7,3],bulkData[10,3],bulkData[13,3])
dataA_adar1b2 <- c(bulkData[1,4],bulkData[4,4],bulkData[7,4],bulkData[10,4],bulkData[13,4])
time <- c(1,2,3,4,5)
pdf('bulk.pdf')
plot(time,dataA_c,type='b',ylim=c(0,900),col="red",cex=1,pch=16,lwd=2,ylab=NULL)
par(new=TRUE)
plot(time,dataA_adar1,type='b',ylim=c(0,900),col="darkblue",cex=1,pch=16,lwd=2,ylab=NULL)
par(new=TRUE)
plot(time,dataA_adarb2,type='b',ylim=c(0,900),col="DarkTurquoise",cex=1,pch=16,lwd=2,ylab=NULL)
par(new=TRUE)
plot(time,dataA_adar1b2,type='b',ylim=c(0,900),col="yellowgreen",cex=1,pch=16,lwd=2,ylab='tumor size')
legend(1,900,c('C','adar1','adarb2','adar12'),col=c('red',"darkblue","DarkTurquoise","yellowgreen"),text.col=c("red","darkblue","DarkTurquoise","yellowgreen") )
dev.off()

A6 <- A_bulkData[1:5,]
A9 <- A_bulkData[6:10,]
A13 <- A_bulkData[11:15,]
A15 <- A_bulkData[16:20,]
A27 <- A_bulkData[21:25,]

t.test(A6['C_bulk'],A6['adar1_bulk'], var.equal=TRUE)
t.test(A6['C_bulk'],A6['adarb2_bulk'], var.equal=TRUE)
t.test(A6['C_bulk'],A6['adar12_bulk'], var.equal=TRUE)

t.test(A9['C_bulk'],A9['adar1_bulk'], var.equal=TRUE)
t.test(A9['C_bulk'],A9['adarb2_bulk'], var.equal=TRUE)
t.test(A9['C_bulk'],A9['adar12_bulk'], var.equal=TRUE)

t.test(A13['C_bulk'],A13['adar1_bulk'], var.equal=TRUE)
t.test(A13['C_bulk'],A13['adarb2_bulk'], var.equal=TRUE)
t.test(A13['C_bulk'],A13['adar12_bulk'], var.equal=TRUE)

t.test(A15['C_bulk'],A15['adar1_bulk'], var.equal=TRUE)
t.test(A15['C_bulk'],A15['adarb2_bulk'], var.equal=TRUE)
t.test(A15['C_bulk'],A15['adar12_bulk'], var.equal=TRUE)

t.test(A27['C_bulk'],A27['adar1_bulk'], var.equal=TRUE)
t.test(A27['C_bulk'],A27['adarb2_bulk'], var.equal=TRUE)
t.test(A27['C_bulk'],A27['adar12_bulk'], var.equal=TRUE)

A_C <- tumor_weight[1:5,1]
A_adar1 <- tumor_weight[1:5,2]
A_adarb2 <- tumor_weight[1:5,3]
A_adar1b2 <- tumor_weight[1:5,4]







