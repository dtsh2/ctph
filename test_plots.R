require(graphics); require(grDevices)

#setwd("C:/Users/David Hayman/Dropbox/marsden/ORION__sample_results")

library(ggplot2)
library(plyr)
library(reshape2)

data<-read.csv("06032015_run_mfp_scan_v5.csv",header=T,row.names=1)
colnames(data) <- c(#"UniqueID",
  "human_3_S3_R1",
  "human_3-Pur_S17_R1",
  "human_3_S3_R2",
  "human_3-Pur_S17_R2",
  "human_2_S2_R1",
  "human_2_S2_R2",
  "human_2-Pur_S16_R1",
  "human_2-Pur_S16_R2",
  "Undetermined_S0_R2",
  "Undetermined_S0_R1",
  "gorilla_4_S4_R1",
  "gorilla_4-Pur_S18_R1",
  "gorilla_4_S4_R2",
  "gorilla_4-Pur_S18_R2",
  "cattle_60-Pur_S26_R1",
  "cattle_60-Pur_S26_R2",
  "cattle_60_S12_R1",
  "cattle_60_S12_R2",
  "cattle_56_S11_R1",
  "cattle_56_S11_R2",
  "cattle_56-Pur_S25_R1",
  "cattle_56-Pur_S25_R2",
  "cattle_62-Pur_S27_R2",
  "cattle_62-Pur_S27_R1",
  "cattle_62_S13_R1",
  "cattle_62_S13_R2",
  "human_1_S1_R2",
  "human_1_S1_R1",
  "human_1-Pur_S15_R1",
  "human_1-Pur_S15_R2",
  "gorilla_44_S8_R1",
  "gorilla_44_S8_R2",
  "gorilla_44-Pur_S22_R1",
  "gorilla_44-Pur_S22_R2",
  "gorilla_24-Pur_S20_R1",
  "gorilla_24_S6_R1",
  "gorilla_24_S6_R2",
  "gorilla_24-Pur_S20_R2",
  "gorilla_14-Pur_S19_R2",
  "gorilla_14-Pur_S19_R1",
  "gorilla_14_S5_R2",
  "gorilla_14_S5_R1",
  "gorilla_34-Pur_S21_R2",
  "gorilla_34-Pur_S21_R1",
  "gorilla_34_S7_R2",
  "gorilla_34_S7_R1",
  "cattle_52_S10_R2",
  "cattle_52_S10_R1",
  "cattle_52-Pur_S24_R2",
  "cattle_52-Pur_S24_R1",
  "cattle_48_S9_R2",
  "cattle_48_S9_R1",	
  "cattle_48-Pur_S23_R2",	
  "cattle_48-Pur_S23_R1",
  "Ctl_S14_R1",
  "Ctl_S14_R2",
  "Ctl-Pur_S28_R1",
  "Ctl-Pur_S28_R2")

data<-data.matrix(data)

## to remove the contamination
df1<-data[!(rowSums(data[,55:58])>0),]
heatmap(df1,Colv=NA,labRow=NA)

test<-df1[,1:54]
# heatmap(test,Colv=NA,labRow=NA,margins=c(12,7),cexCol=0.5)
colnames = c( "human_3",
              "human_2",
              "gorilla_4",
              "cattle_60",
              "cattle_56",
              "cattle_62",
              "human_1",
              "gorilla_44",
              "gorilla_24",
              "gorilla_14",
              "gorilla_34",
              "cattle_52",
              "cattle_48")

# calculate the means
means = lapply(colnames, function(name) { apply(test[,grep(name, colnames(test))], 1, mean) })

# build the result
result = do.call(cbind, means)
#result = as.data.frame(t(result))
colnames(result) = colnames

m <- result[,order(colnames(result))]

rowRm<-which(m[1:6,]==0)

m2<-m[-c(rowRm),]

heatmap(m,Colv=NA,
  labRow=NA,col=terrain.colors(20, alpha = 1),
         margins=c(12,7),
         ColSideColors = c(rep("grey80",5),
                           rep("grey45",5),
                           rep("grey10",3)),
         cexCol=1.8)

heatmap(m2, labRow=NA,col=terrain.colors(20, alpha = 1),
        margins=c(12,7),
        ColSideColors = c(rep("grey45",5),
                          rep("grey80",5),
                          rep("grey10",3)),
        cexCol=1.8)

res<-rownames(test)
 # write.csv(res,"res.csv")

virusresults<-test[grep("virus",res),] # returns those with 'virus' in

virusres<-as.data.frame(virusresults)
 virusres$virus<-rownames(virusres)

virusres<-virusres[-c(1,6),]

vres<-melt(virusres)

virusres[virusres == 0]<-NA
virusres <- virusres[,colSums(is.na(virusres))<nrow(virusres)]
vres<-melt(virusres)

names(vres)<-c("virus","host","value")

colnames = c( "human_3",
              "human_2",
              "gorilla_4",
              "cattle_60",
              "cattle_56",
              "cattle_62",
              "human_1",
              "gorilla_44",
              "gorilla_24",
              "gorilla_14",
              "gorilla_34",
              "cattle_52",
              "cattle_48")

# calculate the means
means = lapply(colnames, function(name) { apply(virusresults[,grep(name, colnames(virusresults))], 1, mean) })

# build the result
result = do.call(cbind, means)
result = as.data.frame(t(result))
rownames(result) = colnames

### KEEP TO REMOVE HAs IF NEEDED
result[result == 0]<-NA
result <- result[rowSums(is.na(result))<ncol(result),]
result

colnames(result)<-c("Microvirus","Picobirnavirus",
                    "D - picobirnavirus","Hepatitis C virus",
                    "Lassa virus","Hyposoter fugitivus ichnovirus")
result<-result[,2:5]

vs_r<-melt(t(result))
colnames(vs_r)<-c("virus","host","reads")

p <- qplot(x=virus, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="dodge",ylab="reads")

p+theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")

p <- qplot(x=virus, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="stack",ylab="reads",asp=1)

p+theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")

ggplot(data=vs_r, aes(x=host,y=reads, fill=host)) +
    labs(y="Mean reads") +
    geom_bar(position='dodge',stat="identity")+
    facet_wrap( ~virus , nrow=1)+
   theme( axis.text.x = element_text(angle = 90,colour = "grey10"),
     text = element_text(size=20))

  colnames(vs_r)[2]<-c("host_id")
  ggplot(data=vs_r, aes(x=host_id,y=reads, fill=host_id)) +
    labs(y="mean reads") +
    geom_bar(position='dodge',stat="identity")+
    facet_wrap( ~virus , nrow=1)+
    theme( axis.text.x = element_blank(),
           axis.title.x = element_blank(),
           text = element_text(size=20))

##############################################
## salmonella
res<-rownames(test)
# write.csv(res,"res.csv")
##
salresults<-test[grep("Salmonella",res),] # returns those with 'Salmonella' in

salres<-as.data.frame(salresults)
salres$sal<-rownames(salres)

salres<-salres[-c(1,6),]

sres<-melt(salres)


salres[salres == 0]<-NA
salres <- salres[,colSums(is.na(salres))<nrow(salres)]
#salres
sres<-melt(salres)

##

names(sres)<-c("sal","host","value")


# calculate the means
means = lapply(colnames, function(name) { apply(salresults[,grep(name, colnames(salresults))], 1, mean) })

# build the result
result = do.call(cbind, means)
result = as.data.frame(t(result))
rownames(result) = colnames

### KEEP TO REMOVE HAs IF NEEDED
result[result == 0]<-NA
result <- result[rowSums(is.na(result))<ncol(result),]

## CHANGE THESE #####
names(result)


colnames(result)<-c("S Enteritidis str. 21180",
"S Enteritidis str. 1753",
"S enterica subsp. enterica",        
"S Agona str. SL483",   
"S Enteritidis str. 84644",
"S enterica subsp. arizonae 62")
#result<-result[,2:5]

vs_r<-melt(t(result))
colnames(vs_r)<-c("Salmonella","host","reads")


p <- qplot(x=Salmonella, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="dodge",ylab="reads")

p+theme(text = element_text(size=15),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")


p <- qplot(x=Salmonella, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="stack",ylab="reads",asp=1)

p+theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")


ggplot(data=vs_r, aes(x=host,y=reads, fill=host)) +
  labs(y="Mean reads") +
  geom_bar(position='dodge',stat="identity")+
  facet_wrap( ~Salmonella , nrow=1)+
  theme( axis.text.x = element_text(angle = 90,colour = "grey10"),
         text = element_text(size=20))

colnames(vs_r)[2]<-c("host_id")
ggplot(data=vs_r, aes(x=host_id,y=reads, fill=host_id)) +
  labs(y="mean reads") +
  geom_bar(position='dodge',stat="identity")+
  facet_wrap( ~Salmonella , nrow=1)+
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         text = element_text(size=20),
         strip.text = element_text(size=15,angle = 90))

##############################################
## influenza
res<-rownames(test)
# write.csv(res,"res.csv")
##
fluresults<-test[grep("influenza",res),] # returns those with 'Salmonella' in

flures<-as.data.frame(fluresults)
flures$flu<-rownames(flures)

#flures<-flures[-c(1,6),]

fres<-melt(flures)

flures[flures == 0]<-NA
flures <- flures[,colSums(is.na(flures))<nrow(flures)]
flures
fres<-melt(flures)


names(fres)<-c("flu","host","value")

# calculate the means
means = lapply(colnames, function(name) { apply(fluresults[,grep(name, colnames(fluresults))], 1, mean) })

# build the result
result = do.call(cbind, means)
result = as.data.frame(t(result))
rownames(result) = colnames

### KEEP TO REMOVE HAs IF NEEDED
result[result == 0]<-NA
result <- result[rowSums(is.na(result))<ncol(result),]
# result

## CHANGE THESE #####
names(result)


colnames(result)<-c("H. parainfluenzae T3T1",
                    "H. influenzae")
#result<-result[,2:5]

vs_r<-melt(t(result))
colnames(vs_r)<-c("Haemophilus","host","reads")


p <- qplot(x=Haemophilus, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="dodge",ylab="reads")

p+theme(text = element_text(size=15),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")


p <- qplot(x=Haemophilus, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="stack",ylab="reads",asp=1)

p+theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")


ggplot(data=vs_r, aes(x=host,y=reads, fill=host)) +
  labs(y="Mean reads") +
  geom_bar(position='dodge',stat="identity")+
  facet_wrap( ~Haemophilus , nrow=1)+
  theme( axis.text.x = element_text(angle = 90,colour = "grey10"),
         text = element_text(size=20))

colnames(vs_r)[2]<-c("host_id")
ggplot(data=vs_r, aes(x=host_id,y=reads, fill=host_id)) +
  labs(y="mean reads") +
  geom_bar(position='dodge',stat="identity")+
  facet_wrap( ~Haemophilus , nrow=1)+
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         text = element_text(size=20),
         strip.text = element_text(size=15))

#######################
## plasmodium
res<-rownames(test)
# write.csv(res,"res.csv")
##
plasresults<-test[grep("Plasmodium",res),] # returns those with 'Salmonella' in

plasres<-as.data.frame(plasresults)
plasres$plas<-rownames(plasres)

#plasres<-plasres[-c(1,6),]

pres<-melt(plasres)

plasres[plasres == 0]<-NA
plasres <- plasres[,colSums(is.na(plasres))<nrow(plasres)]
pres<-melt(plasres)

names(pres)<-c("plas","host","value")


# calculate the means
means = lapply(colnames, function(name) { apply(plasresults[,grep(name, colnames(plasresults))], 1, mean) })

# build the result
result = do.call(cbind, means)
result = as.data.frame(t(result))
rownames(result) = colnames

### KEEP TO REMOVE HAs IF NEEDED
result[result == 0]<-NA
result <- result[rowSums(is.na(result))<ncol(result),]
#result

## CHANGE THESE #####
names(result)


colnames(result)<-c("P. reichenowi",                 
"P. knowlesi strain H",
"P. falciparum NF135/5.C10",
"P. falciparum",
"P. falciparum Palo Alto/Uganda",
"P. falciparum Tanzania")
#result<-result[,2:5]

vs_r<-melt(t(result))
colnames(vs_r)<-c("Plasmodium","host","reads")


p <- qplot(x=Plasmodium, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="dodge",ylab="reads")

p+theme(text = element_text(size=15),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")


p <- qplot(x=Plasmodium, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="stack",ylab="reads",asp=1)

p+theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")


ggplot(data=vs_r, aes(x=host,y=reads, fill=host)) +
  labs(y="Mean reads") +
  geom_bar(position='dodge',stat="identity")+
  facet_wrap( ~Plasmodium , nrow=1)+
  theme( axis.text.x = element_text(angle = 90,colour = "grey10"),
         text = element_text(size=20))

colnames(vs_r)[2]<-c("host_id")
ggplot(data=vs_r, aes(x=host_id,y=reads, fill=host_id)) +
  labs(y="mean reads") +
  geom_bar(position='dodge',stat="identity")+
  facet_wrap( ~Plasmodium , nrow=1)+
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         text = element_text(size=20),
         strip.text = element_text(size=15,angle=90))
##
####### Trypanosoma

res<-rownames(test)
# write.csv(res,"res.csv")
##
trypresults<-test[grep("Trypanosoma",res),] # returns those with 'Salmonella' in

trypres<-as.data.frame(trypresults)
trypres$tryp<-rownames(trypres)

#trypres<-trypres[-c(1,6),]

tryres<-melt(trypres)
# 

trypres[trypres == 0]<-NA
trypres <- trypres[,colSums(is.na(trypres))<nrow(trypres)]
# trypres
tryres<-melt(trypres)

tryres$trypanosome<-rep(c("T. congolense",
                    "T. vivax",
                    "T. brucei gambiense"),8)

p <- qplot(x=trypanosome, y=value, fill=variable,
           data=tryres, geom="bar", stat="identity",
           position="dodge",ylab="reads")

p+theme(text = element_text(size=18),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")

##

names(tryres)<-c("tryp","host","value")

# calculate the means
means = lapply(colnames, function(name) { apply(trypresults[,grep(name, colnames(trypresults))], 1, mean) })

# build the result
result = do.call(cbind, means)
result = as.data.frame(t(result))
rownames(result) = colnames

### KEEP TO REMOVE HAs IF NEEDED
result[result == 0]<-NA
result <- result[rowSums(is.na(result))<ncol(result),]
#result

## CHANGE THESE #####
names(result)


colnames(result)<-c("T. congolense",
"T. vivax",
"T. brucei gambiense")
#result<-result[,2:5]

vs_r<-melt(t(result))
colnames(vs_r)<-c("Trypanosoma","host","reads")

p <- qplot(x=Trypanosoma, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="dodge",ylab="reads")

p+theme(text = element_text(size=15),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")

p <- qplot(x=Trypanosoma, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="stack",ylab="reads",asp=1)

p+theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")

ggplot(data=vs_r, aes(x=host,y=reads, fill=host)) +
  labs(y="Mean reads") +
  geom_bar(position='dodge',stat="identity")+
  facet_wrap( ~Trypanosoma , nrow=1)+
  theme( axis.text.x = element_text(angle = 90,colour = "grey10"),
         text = element_text(size=20))

colnames(vs_r)[2]<-c("host_id")
ggplot(data=vs_r, aes(x=host_id,y=reads, fill=host_id)) +
  labs(y="mean reads") +
  geom_bar(position='dodge',stat="identity")+
  facet_wrap( ~Trypanosoma , nrow=1)+
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         text = element_text(size=20),
         strip.text = element_text(size=15,angle=90))

##
####### Coxiella

res<-rownames(test)
# write.csv(res,"res.csv")
##
coxresults<-test[grep("Coxiella",res),] # returns those with 'Salmonella' in
coxresults<-as.data.frame(coxresults)
colnames(coxresults)<-"360115:Coxiella burnetii RSA 331"
coxresults<-as.data.frame(t(coxresults))

 coxres<-melt(coxresults)
 coxres$cox<-rownames(coxresults)

coxres[coxres == 0]<-NA
coxres <- coxres[!(is.na(coxres$value)),]
names(coxres)<-c("host","value","cox")

# calculate the means
means = lapply(colnames, function(name) { apply(coxresults[,grep(name, colnames(coxresults))], 1, mean) })
means<-as.data.frame(means)
# build the result
result = do.call(cbind, means)
result = as.data.frame(t(result))
rownames(result) = colnames

##########################
colnames(result)<-c("C. burnetii RSA")
result$host<-rownames(result)
result$Coxiella<-'C. burnetti'

result<-result[c(4:7,13),]

vs_r<-result
colnames(vs_r)<-c("reads","host","Coxiella")

p <- qplot(x=Coxiella, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="dodge",ylab="reads")

p+theme(text = element_text(size=15),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")

ggplot(data=vs_r, aes(x=host,y=reads, fill=host)) +
  labs(y="Mean reads") +
  geom_bar(position='dodge',stat="identity")+
  facet_wrap( ~Coxiella , nrow=1)+
  theme( axis.text.x = element_text(angle = 90,colour = "grey10"),
         text = element_text(size=20))

colnames(vs_r)[2]<-c("host_id")
ggplot(data=vs_r, aes(x=host_id,y=reads, fill=host_id)) +
  labs(y="mean reads") +
  geom_bar(position='dodge',stat="identity")+
  facet_wrap( ~Coxiella , nrow=1)+
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         text = element_text(size=20),
         strip.text = element_text(size=15))

####### Treponema

res<-rownames(test)
# write.csv(res,"res.csv")
##
trepresults<-test[grep("Treponema",res),] # returns those with 'Salmonella' in

trepres<-as.data.frame(trepresults)
trepres$trep<-rownames(trepres)

#trepres<-trepres[-c(1,6),]

treres<-melt(trepres)
# 

trepres[trepres == 0]<-NA
trepres <- trepres[,colSums(is.na(trepres))<nrow(trepres)]
# trepres
treres<-melt(trepres)

treres$treponema<-rep(c("T. primitia ZAS-2",
"T. azotonutricium ZAS-9",
"T. caldaria DSM",
"T. sp. OMZ",
"T. pedis str. T",
"T. denticola"),20)

p <- qplot(x=treponema, y=value, fill=variable,
           data=treres, geom="bar", stat="identity",
           position="dodge",ylab="reads")

p+theme(text = element_text(size=18),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")

##

names(treres)<-c("trep","host","value")

# calculate the means
means = lapply(colnames, function(name) { apply(trepresults[,grep(name, colnames(trepresults))], 1, mean) })

# build the result
result = do.call(cbind, means)
result = as.data.frame(t(result))
rownames(result) = colnames

### KEEP TO REMOVE HAs IF NEEDED
result[result == 0]<-NA
result <- result[rowSums(is.na(result))<ncol(result),]
#result

## CHANGE THESE #####
names(result)


colnames(result)<-c("T. primitia ZAS-2",
"T. azotonutricium ZAS-9",
"T. caldaria DSM",
"T. sp. OMZ",
"T. pedis str. T",
"T. denticola")
#result<-result[,2:5]

vs_r<-melt(t(result))
colnames(vs_r)<-c("Treponema","host","reads")

p <- qplot(x=Treponema, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="dodge",ylab="reads")

p+theme(text = element_text(size=15),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")

p <- qplot(x=Treponema, y=reads, fill=host,
           data=vs_r, geom="bar", stat="identity",
           position="stack",ylab="reads",asp=1)

p+theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, vjust=1,colour="grey20"))+
  scale_fill_hue(c=90, l=55)+ scale_fill_brewer(palette="Set1")

ggplot(data=vs_r, aes(x=host,y=reads, fill=host)) +
  labs(y="Mean reads") +
  geom_bar(position='dodge',stat="identity")+
  facet_wrap( ~Treponema , nrow=1)+
  theme( axis.text.x = element_text(angle = 90,colour = "grey10"),
         text = element_text(size=20))

colnames(vs_r)[2]<-c("host_id")
ggplot(data=vs_r, aes(x=host_id,y=reads, fill=host_id)) +
  labs(y="mean reads") +
  geom_bar(position='dodge',stat="identity")+
  facet_wrap( ~Treponema , nrow=1)+
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         text = element_text(size=20),
         strip.text = element_text(size=15,angle=90))
