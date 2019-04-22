rm(list=ls())

# --------------------------------------------------------------------
# Code de Matthieu pour fusionner les fragments proches 
# --------------------------------------------------------------------

par(mfrow=c(1,1))
data = read.table("8F0EB0CW11N-Alignment.txt", sep='\t', header = TRUE)
x = TRUE
while (x==TRUE){
  x = FALSE
  for (i in 1:length(data$q._start)) {
    if (data$alignment_length[i]>=3000) {
      for (j in 1:length(data$q._start)) {
        if (i!=j){
          prem = abs(data$q._end[i]-data$q._start[j])
          deux = abs(data$s._end[i]-data$s._start[j])
          troi = abs(data$q._end[j]-data$q._start[i])
          quat = abs(data$s._end[j]-data$s._start[i])
          if (prem<=15000 && deux<=15000 && ((data$s._start[i]-data$s._end[i]>0 && data$s._start[j]-data$s._end[j]>0) || (data$s._start[i]-data$s._end[i]<0 && data$s._start[j]-data$s._end[j]<0))){
            data$q._end[i]=data$q._end[j]
            data$s._end[i]=data$s._end[j]
            data$alignment_length[i]=max(abs(data$q._end[i]-data$q._start[i]),abs(data$s._end[i]-data$s._start[i]))
            data[j,]=0
            x = TRUE
          }
          else if (troi<=15000 && quat<=15000 && ((data$s._start[i]-data$s._end[i]>0 && data$s._start[j]-data$s._end[j]>0) || (data$s._start[i]-data$s._end[i]<0 && data$s._start[j]-data$s._end[j]<0))){
            data$q._start[i]=data$q._start[j]
            data$s._start[i]=data$s._start[j]
            data$alignment_length[i]=max(abs(data$q._end[i]-data$q._start[i]),abs(data$s._end[i]-data$s._start[i]))
            data[j,]=0
            x = TRUE
          }
        }
      }
    }
  }
}
# --------------------------------------------------------------------
# Selection des fragments longs
# --------------------------------------------------------------------
orthologues=subset(data,alignment_length>=50000)
autres=subset(data,alignment_length<50000)

# --------------------------------------------------------------------
# Identification des inversions
# --------------------------------------------------------------------
inversion = orthologues$q._start > orthologues$q._end | orthologues$s._start > orthologues$s._end
couleur = ifelse(inversion == TRUE, "red" , "green" )
orthologues=cbind(orthologues,couleur)

# --------------------------------------------------------------------
# Dot plot filtré
# --------------------------------------------------------------------
plot(c(orthologues$q._start[1],orthologues$q._end[1]),c(orthologues$s._start[1],orthologues$s._end[1]),xlim=c(0,4500000),
     ylim=c(0,4500000),'l',xlab="Génome médiéval", ylab="Yersinia pestis CO92",
     main="Dotplot filtré",col="white")

for (l in 1:length(orthologues$q._start)) {
  print(orthologues$couleur[l])
  segments(orthologues$q._start[l],orthologues$s._start[l],orthologues$q._end[l],orthologues$s._end[l], col = orthologues$couleur[l])
}

# --------------------------------------------------------------------
# Stockage dans un fichier à part
# --------------------------------------------------------------------
write.table(orthologues,"8F0EB0CW11N-fragments.csv",row.names=TRUE,sep=",",dec=".", na=" ")

sel = data > 5000
table2 = table1[sel,c("alignment_length","q._start","q._end","s._start","s._end") ]
table3 = table1[!sel,]
dim(table2) # taille filtrée => orthologues
dim(table3) # non orthologues
head(table2)
tail(table2)

mean(orthologues$X._identity)  
mean(autres$X._identity)
mean(data$X._identity)

# Trier le tableau selon ordre croissant q.start
tri_med = orthologues[,c("alignment_length","q._start","q._end","s._start","s._end") ]
tri_med <- tri_med[order(tri_med$q._start),]


# --------------------------------------------------------------------
# Creation d'une colonne EDGE entre les gènes pour les 2 génomes
# --------------------------------------------------------------------
v_num = vector()
edges_med = vector()
for (i in 1:nrow(tri_med)) {
  if (tri_med$inversion[i]  == TRUE) {
    v_num[i] =paste("med",i,"inv",sep="_") }
  else {v_num[i] =paste("med",i,sep="_") }
  
  if (i == nrow(tri_med)) {
    edges_med[i] =paste(tri_med$q._end[i], '\t', tri_med$q._start[1])
  }
  
  else{
  if (tri_med$inversion[i]  == FALSE & tri_med$inversion[i+1]  == FALSE) {
    edges_med[i] =paste(tri_med$q._end[i], '\t', tri_med$q._start[i+1]) } 
  
  if (tri_med$inversion[i]  == TRUE & tri_med$inversion[i+1]  == TRUE) {
    edges_med[i] =paste(tri_med$q._start[i], '\t', tri_med$q._end[i+1]) }

  if (tri_med$inversion[i]  == TRUE & tri_med$inversion[i+1]  == FALSE ) {
    edges_med[i] =paste(tri_med$q._start[i], '\t', tri_med$q._start[i+1])
    print("TRUE FALSE")
    print(i)
    print(edges_med[i])}
  
  if (tri_med$inversion[i]  == FALSE & tri_med$inversion[i+1]  == TRUE ) {
    edges_med[i] =paste(tri_med$q._end[i], '\t', tri_med$q._end[i+1])
    print("FALSE TRUE")
    print(i)
    print(edges_med[i])}
  }
  }


tri_med = cbind(tri_med, v_num)
tri_med = cbind(tri_med, edges_med)

# Trier le tableau selon ordre croissant s.start
tri_CO92 <- tri_med[order(tri_med$s._start),]

v_num2 = vector()
edges_CO92 = vector() 

for (i in 1:nrow(tri_CO92)) {
  if (tri_CO92$inversion[i]  == TRUE) {
    v_num2[i] =paste("CO92",i,"inv",sep="_") }
  else {v_num2[i] =paste("CO92",i,sep="_")}

  if (i == nrow(tri_CO92)) {
    edges_CO92[i] =paste(tri_CO92$q._end[i], '\t', tri_CO92$q._start[1])
  }
 else{ 
  if (tri_CO92$inversion[i]  == FALSE & tri_CO92$inversion[i+1]  == FALSE) {
    edges_CO92[i] =paste(tri_CO92$q._end[i], '\t', tri_CO92$q._start[i+1]) } 

  if (tri_CO92$inversion[i]  == TRUE & tri_CO92$inversion[i+1]  == TRUE) {
   edges_CO92[i] =paste(tri_CO92$q._start[i], '\t', tri_CO92$q._end[i+1]) }

  if (tri_CO92$inversion[i]  == TRUE & tri_CO92$inversion[i+1]  == FALSE ) {
    edges_CO92[i] =paste(tri_CO92$q._start[i], '\t', tri_CO92$q._start[i+1])
}

  if (tri_CO92$inversion[i]  == FALSE & tri_CO92$inversion[i+1]  == TRUE ) {
   edges_CO92[i] =paste(tri_CO92$q._end[i], '\t', tri_CO92$q._end[i+1])
 }
}
}

final = cbind(tri_CO92, v_num2)
final = cbind(tri_CO92, edges_CO92)

head(tri_med)
head(tri_CO92)
head(final)

# Tri avant écriture
final = final[order(final$q._start),]
write.table(final, "final.csv", row.names=FALSE, sep=",",dec=".", na=" ")
 