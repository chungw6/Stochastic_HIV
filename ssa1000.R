
library(dplyr)

dat1 = read.table("/Users/cwjia/Documents/Virus Work/c_ssa1000_2.txt")
colnames(dat1) = c("RunNo","First","Compartment","EndTime","Time")
sapply(dat1[,1:5], mean)
head(dat1)
summary(dat1)

Y1e = dat1[dat1$Compartment==1,]
Y1f = dat1[dat1$Compartment==2,]

cols <- c("First","EndTime","Time")

dim(Y1e)
dim(Y1f)

summary(Y1e[cols])
summary(Y1f[cols])

sapply(Y1e[cols], sd)
sapply(Y1f[cols], sd)
