rm(list=ls())

par(mfrow=c(1,1))
data=read.table("8EZHDCFD11N-Alignment.txt",dec=".",h=F)
data = read.table("8F0EB0CW11N-Alignment.txt", sep='\t', header = TRUE)
x = TRUE
while (x==TRUE){
  x = FALSE
  for (i in 1:length(data$V7)) {
    if (data$V4[i]>=3000) {
      for (j in 1:length(data$V7)) {
        if (i!=j){
          prem = abs(data$V8[i]-data$V7[j])
          deux = abs(data$V10[i]-data$V9[j])
          troi = abs(data$V8[j]-data$V7[i])
          quat = abs(data$V10[j]-data$V9[i])
          if (prem<=15000 && deux<=15000 && ((data$V9[i]-data$V10[i]>0 && data$V9[j]-data$V10[j]>0) || (data$V9[i]-data$V10[i]<0 && data$V9[j]-data$V10[j]<0))){
            data$V8[i]=data$V8[j]
            data$V10[i]=data$V10[j]
            data$V4[i]=max(abs(data$V8[i]-data$V7[i]),abs(data$V10[i]-data$V9[i]))
            data[j,]=0
            x = TRUE
          }
          else if (troi<=15000 && quat<=15000 && ((data$V9[i]-data$V10[i]>0 && data$V9[j]-data$V10[j]>0) || (data$V9[i]-data$V10[i]<0 && data$V9[j]-data$V10[j]<0))){
            data$V7[i]=data$V7[j]
            data$V9[i]=data$V9[j]
            data$V4[i]=max(abs(data$V8[i]-data$V7[i]),abs(data$V10[i]-data$V9[i]))
            data[j,]=0
            x = TRUE
          }
        }
      }
    }
  }
}

inversion = data$V9 > data$V10
table(inversion)
couleur = ifelse(inversion == TRUE, "Red" , "Green" )
table(couleur)

plot(c(data$V7[1],data$V8[1]),c(data$V9[1],data$V10[1]),xlim=c(0,4500000),
     ylim=c(0,4500000),'l',xlab="Génome médiéval", ylab="Yersinia pestis CO92",
     main="Dotplot filtré")

for (l in 2:length(data$V7)) {
  if (data$V4[l]>=50000){
    segments(data$V7[l],data$V9[l],data$V8[l],data$V10[l],col = couleur[l])
  }
}

newdata=subset(data,V4>=50000)
inversion = newdata$V7 > newdata$V8 | newdata$V9 > newdata$V10
couleur = ifelse(inversion == TRUE, "red" , "green" )
newdata=cbind(newdata,couleur)
plot(c(newdata$V7[1],newdata$V8[1]),c(newdata$V9[1],newdata$V10[1]),xlim=c(0,4500000),
     ylim=c(0,4500000),'l',xlab="Génome médiéval", ylab="Yersinia pestis CO92",
     main="Dotplot filtré",col="white")

for (l in 1:length(newdata$V7)) {
  print(newdata$couleur[l])
  segments(newdata$V7[l],newdata$V9[l],newdata$V8[l],newdata$V10[l], col = newdata$couleur[l])
}
write.table(newdata,"8F0EB0CW11N-fragments.txt",sep="\t",row.names=TRUE)
