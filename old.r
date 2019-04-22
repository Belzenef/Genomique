table1 = read.table("8F0EB0CW11N-Alignment.txt", sep='\t', header = TRUE)
head(table1)
names(table1)
dim(table1)

# query_acc.ver = génome médiéval 1348
# subject.acc.ver = Yersinia pestis CO92 NC_003143.1

sel = table1$alignment_length > 3000
table2 = table1[sel,c("alignment_length","q._start","q._end","s._start","s._end") ]
table3 = table1[!sel,]
dim(table2) # taille filtrée => orthologues
dim(table3) # non orthologues
head(table2)
tail(table2)


plot(table2[1,"q._start"], table2[1,"s._start"], main="Dotplot filtré",
     xlab="Génome médiéval", ylab="Yersinia pestis CO92",
     xlim=c(0, max(table2$q._end)), ylim=c(0, max(table2$s._end)), col="white")

 for (i in 1:nrow(table2)){
  segments(table2[i,"q._start"], table2[i,"s._start"], table2[i,"q._end"], 
           table2[i,"s._end"], col =i, lwd=2)
}

inversion = table2["s._start"] > table2["s._end"]
table(inversion)
table2$inversion = inversion

couleur = ifelse(table2$inversion == TRUE, "Red" , "Green" )

# Segment vert = non inversés, rouge = inversion  
for (i in 1:nrow(table2)){
  segments(table2[i,"q._start"], table2[i,"s._start"], table2[i,"q._end"], 
           table2[i,"s._end"], col = couleur[i], lwd=2) }

mean(table1$X._identity)  
mean(table1[sel,]$X._identity)
mean(table1[-sel,]$X._identity)

# Trier le tableau selon ordre croissant q.start
tri_med <- table2[order(table2$q._start),]

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
}
}


final = cbind(tri_CO92, v_num2)
final = cbind(tri_CO92, edges_CO92)

head(tri_med)
head(tri_CO92)
head(final)

# Tri avant écriture
final = final[order(final$q._start),]
write.table(final, "final2.csv", row.names=FALSE, sep=",",dec=".", na=" ")
 