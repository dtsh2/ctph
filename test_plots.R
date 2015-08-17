require(graphics); require(grDevices)

setwd("C:/Users/David Hayman/Dropbox/marsden/ORION__sample_results")

# data<-read.csv("06032015_run_mfp_scan_v5.csv",header=T,row.names=1)
# 
# data<-data.matrix(data)
# heatmap(data)
# heatmap(data,na.rm=T)

data<-read.csv("marsden_test.csv",header=T,row.names=1)
# data<-read.csv("marsden_test_1.csv",header=T,row.names=1)
# data<-read.csv("marsden_test_2.csv",header=T,row.names=1)
# data<-read.csv("marsden_test_3.csv",header=T,row.names=1)
data<-data.matrix(data)
# heatmap(data)
# heatmap(data,na.rm=T)

# df2 <- cbind(data[,1]-data[,10],
#              data[,2]-data[,10],
#              data[,3]-data[,10],
#              data[,4]-data[,10],
#              data[,5]-data[,10],
#              data[,6]-data[,10],
#              data[,7]-data[,10],
#              data[,8]-data[,10],
#              data[,9]-data[,10],
#              data[,10]-data[,10])
# df2 <- cbind(data[,1]-data[,4],
#              data[,2]-data[,4],
#              data[,3]-data[,4],
#              data[,4]-data[,4])

df2 <- cbind(data[,1]-data[,14],
             data[,2]-data[,14],
             data[,3]-data[,14],
             data[,4]-data[,14],
             data[,5]-data[,14],
             data[,6]-data[,14],
             data[,7]-data[,14],
             data[,8]-data[,14],
             data[,9]-data[,14],
             data[,10]-data[,14],
             data[,11]-data[,14],
             data[,12]-data[,14],
             data[,13]-data[,14],
             data[,14]-data[,14])

df2[df2<0] <- 0
colnames(df2)<-colnames(data)

#heatmap(df2)

df2[df2 == 0] <- NA

df2 <- df2[,colSums(is.na(df2))<nrow(df2)]
df2 <- df2[!rowSums(is.na(df2)) == ncol(df2),]

df2[is.na(df2)] <- 0
#df2<-df2[,1:13]
summary(df2)

# df2<-data.matrix(df2)

# heatmap(df2)

test<-df2[ order(-df2[,1]), ]
test<-test[-c(1),]
# heatmap(test)
# heatmap(test,Colv=NA,labRow="")

# heatmap(test,Colv=NA,
#         labRow="",#revC=T,
#         margins=c(8,7),
#         ColSideColors = c(rep("lightblue",3),
#                            rep("blue",5),
#                            rep("darkblue",5)),
#         cexCol=1.8)
## # 19695 # hits

heatmap(test,Colv=NA,
        labRow=NA,col=terrain.colors(20, alpha = 1),
        margins=c(12,7),
        ColSideColors = c(rep("lightgrey",3),
                          rep("grey",5),
                          rep("darkgrey",5)),
        cexCol=1.8)

res<-rownames(test)
 # write.csv(res,"res.csv")
virusresults<-test[c(300, # rift valley fever virus - 1 gorilla (4), 3 cattle (1:3) 
1016, # picobirnavirus 1 cow (2) # novel - see map : http://www.hindawi.com/journals/bmri/2014/780752/
1021, # Dromedary picobirnavirus 1 cow (2)
1069#, # Lassa virus 1 human (2)
# 1147 # Hyposoter fugitivus ichnovirus
) ,]
##
rownames(virusresults)<-c("RVFV",#"RVF virus",
                          "PV",#"Picobirnavirus",
                          "PV(D)",#"Picobirnavirus (D)",
                          "LFV")#'Lassa virus')#,
                          #"Ichnovirus (Hyposoter fugitivus)")

virusres<-as.data.frame(virusresults)
virusres$virus<-rownames(virusres)
library(ggplot2)
library(plyr)
library(reshape2)

vres<-melt(virusres)

p <- qplot(x=virus, y=value, fill=variable,
                       data=vres, geom="bar", stat="identity",
                       position="dodge",ylab="reads")

p+theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20")) 


p <- qplot(x=virus, y=value, fill=variable,
           data=vres, geom="bar", stat="identity",
           position="stack",ylab="reads")

p+theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20")) 


virusres[virusres == 0]<-NA
virusres <- virusres[,colSums(is.na(virusres))<nrow(virusres)]
virusres
vres<-melt(virusres)

p <- qplot(x=virus, y=value, fill=variable,
           data=vres, geom="bar", stat="identity",
           position="dodge",ylab="reads")

p+theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")

p <- qplot(x=virus, y=value, fill=variable,
           data=vres, geom="bar", stat="identity",
           position="stack",ylab="reads")

# cbPalette <- c("#999999", "#56B4E9", "#0072B2", "#D55E00", "#CC79A7")

p+theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")


##

names(vres)<-c("virus","host","value")

p <- qplot(x=virus, y=value, fill=host,
           data=vres, geom="bar", stat="identity",
           position="dodge",ylab="reads")
p+ theme(text = element_text(size=20),
         axis.title.x=element_blank(),
      axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")

#  ggplot(data=vres, aes(x=variable,y=value, fill=variable)) + 
#    geom_bar(position='dodge')+
#    facet_wrap( ~virus , nrow=1)+
#   theme(axis.title.x=element_blank())

# vres[order(vres$virus),] 
# 
# ##################
# 
# library(gplots)
# heatmap.2(test,dendrogram="row",labRow=NA,
#          col=heat.colors(20),  cexCol=1,
#          scale="col")
# 
# heatmap.2(test,dendrogram="row",
#           # col=heat.colors(120),
#           cexCol=1,scale="row",labRow=NA,
#           density.info="none")
# 
# heatmap.2(test,dendrogram="row", #col=redgreen(75),
#           scale="row", key=T, keysize=1.5,
#           density.info="none", trace="none",cexCol=0.9, labRow=NA)
# 
# heatmap(test,Colv=NA,
#         labRow="",
#         margins=c(8,7),
#         ColSideColors = c(rep("#80FFFFFF",3),
#                           rep("lightgreen",5),
#                           rep("#FF80FFFF",5)),
#         cexCol=2,
#         col=cm.colors(20))
