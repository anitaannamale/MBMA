library(ggplot2)

# Lecture des arguments en ligne de commande
filename = "../data/mock/perc_especes_mdg.txt"

# Lecture des fichiers
d = read.table(filename, header=T ,sep='\t', stringsAsFactors = F)

new_d=c()
Ech = unique(d$Echantillon)
for(i in Ech){
  sample = d[d$Echantillon == i,]
  sample = sample[order(-sample$Percentage),]
  new_d=rbind(new_d,sample)
}

## species
ggplot(new_d,aes(x=new_d$Echantillon,y=new_d$Percentage))+
geom_bar(aes(fill=new_d$Species, width=0.75),colour="white",stat="identity")+
ylab("Pourcentage") +
xlab("Echantillon") +
scale_y_continuous(labels = scales::percent) +
theme(strip.text=element_text(size=14, face='bold'),
      legend.position="None",
      axis.text.x = element_text(angle = 90, hjust=1),
      axis.title = element_text(size = 15, face ='bold'),
      legend.title = element_text(face = 'bold', size = 12))

#### genes de resistance
ggplot(new_d,aes(x=new_d$Echantillon,y=new_d$Percentage))+
  geom_bar(aes(fill=new_d$Resistance.genes, width=0.75),colour="white",stat="identity")+
  ylab("Relative abundance") +
  xlab("Sample") +
  scale_y_continuous(labels = scales::percent) +
  theme(strip.text=element_text(size=14, face='bold'),
        legend.position="None",
        axis.text.x = element_text(angle = 90, hjust=1),
        axis.title = element_text(size = 15, face ='bold'),
        legend.title = element_text(face = 'bold', size = 12))

write.table(new_d, file="resist_krem_ordonné.txt", row.names=FALSE, quote=FALSE, sep='\t', col.names = colnames(new_d))


## calcul de la divergence de KL

# function
KLtest = function(x1,x2)
{
  res = 0
  x1 = x1/sum(x1)
  x2 = x2/sum(x2)
  for(i in 1:length(x1)){
  if(x1[i] != 0 && x2[i] != 0)
  {
    res = res + x1[i]*log(x1[i]/x2[i])
    
  }
  }
  return(res)
}

#value = KLtest(percentage1, percentage2)