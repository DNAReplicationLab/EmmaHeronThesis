##################
# IN UNIX
# replace #h with nothing in the specified cmap and xmap files so that R can read the column headings
# sed 's/#h //g' autoNoise0_q.cmap > autoNoise0_h.cmap
# sed 's/#h //g' autoNoise0.xmap > autoNoise0_h.xmap
#################

# assign the specified cmap and xmap files to the variables PATH_c and PATH_x
PATH_c=("~/path/to/cmap/autoNoise0_h.cmap")
PATH_x=("~/path/to/xmap/autoNoise0_h.xmap")

# specify the cmap variable (data frame) as PATH_c
cmap=data.frame(read.table(PATH_c, header=TRUE, comment.char="#"))

# specify the xmap variable as PATH_x
# had to put a new line ('enter') at the end of the xmap file to complete the final line
xmap=read.table(PATH_x, sep="\t", header=TRUE)

# specify x_cols variable as the vector given (= the columns in xmap numbered 1-8, and 11)
x_cols=c(c(colnames(xmap)[c(seq(from=1,to=8),11)])) 
# specify c_cols variable as the vector given (= the columns in cmap numbered 2, 5, 6, and 10)
c_cols=c(colnames(cmap)[c(2,5,6,10)])

# create a new data frame called "COMPILED" which combines data from the cmap and xmap based on the CMapID/ID
COMPILED=data.frame()
for(ID in unique(cmap$CMapId)){
  print(ID)
  nROW=length(which(cmap$CMapId == ID))
  compiled=cmap[which(cmap$CMapId == ID), c_cols]
  
  xmap_rep=data.frame( apply(xmap[which(xmap$QryContigID==ID), x_cols], 
                             2, FUN=function(x) {rep(x, times=nROW)}))
  
  compiled=data.frame(cbind(compiled, xmap_rep))

  COMPILED=rbind(COMPILED,compiled)
  rm(compiled)
}

#######################
##### FILTERING #######
#######################
# filter out the final row of each molecule which is label channel 0
COMPILED=COMPILED[which(COMPILED$LabelChannel %in% c(1, 2)),]

# change the positional columns to be numerics
positions=c("Position", "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos")
COMPILED[,positions]=apply(COMPILED[,positions], 2, FUN=function(x) {as.numeric(as.character(x))})

# round to nearest base for positional columns
COMPILED[,positions]=apply(COMPILED[,positions], 2, FUN=function(x) {round(x,digits=0)})

# split into + and - orientation molecules
COMPILED_plus=COMPILED[which(COMPILED$Orientation =="+"),]
COMPILED_minus=COMPILED[which(COMPILED$Orientation =="-"),]
write.table(COMPILED_plus,file="cmap_xmap_plus.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(COMPILED_minus,file="cmap_xmap_minus.txt", sep="\t", quote=FALSE, row.names=FALSE)

######################################################
### plus strand chromosomal coordinate calculation ###
######################################################
# calculate distance for each label from the start label
COMPILED_plus$LabelDistance = rep(NA, nrow(COMPILED_plus))
COMPILED_plus$LabelDistance = COMPILED_plus$Position - COMPILED_plus$QryStartPos
COMPILED_plus=COMPILED_plus[which(COMPILED_plus$LabelDistance >= 0),]

# calculate chromosomal coordinate of each label
COMPILED_plus$ChrCoordinate = rep(NA, nrow(COMPILED_plus))
COMPILED_plus$ChrCoordinate = COMPILED_plus$LabelDistance + COMPILED_plus$RefStartPos
write.table(COMPILED_plus,file="compiled_plus.txt", sep="\t", quote=FALSE, row.names=FALSE)

#######################################################
### minus strand chromosomal coordinate calculation ###
#######################################################
# calculate distance for each label from the start label
COMPILED_minus$LabelDistance = rep(NA, nrow(COMPILED_minus))
COMPILED_minus$LabelDistance = COMPILED_minus$Position - COMPILED_minus$QryStartPos

# multiply LabelDistance values by -1
COMPILED_minus$LabelDistanceInvert = rep(NA, nrow(COMPILED_minus))
COMPILED_minus$LabelDistanceInvert = COMPILED_minus$LabelDistance * -1
COMPILED_minus=COMPILED_minus[which(COMPILED_minus$LabelDistanceInvert >= 0),]

# calculate chromosomal coordinate of each label
COMPILED_minus$ChrCoordinate = rep(NA, nrow(COMPILED_minus))
COMPILED_minus$ChrCoordinate = COMPILED_minus$LabelDistanceInvert + COMPILED_minus$RefStartPos
write.table(COMPILED_minus,file="compiled_minus", sep="\t", quote=FALSE, row.names=FALSE)