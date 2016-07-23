# # Lecture des arguments en ligne de commande

library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type = "character", default=NULL, 
              help="count matrice file name", metavar="character"),
  make_option(c("-r", "--reference_db"), type="character", default=NULL,
              help="dabases genes informations", metavar = "character"),
  make_option(c("-l","--label"), type = "character", default = "MBMA",
              help = "label for the analysis [default= %default]", metavar = "character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least two argument must be supplied (input file & reference file)", call.=FALSE)
}

# Lecture des fichiers
matrix = read.table(opt$file, header=T ,sep='\t', stringsAsFactors = F)
ref=read.table(opt$reference_db,header=T, sep='\t')

# Récupération du taxID des souches
TaxID = NULL
for(id in matrix$ref_genome) {
  TaxID = c(TaxID, unlist(strsplit(id, "[.]"))[1])
}

# Attribution des noms de chaque espèces au TaxID
name = NULL
for(id in TaxID){
  species = ref[which(ref$TaxID == id),]$Species
  name=c(name, as.character(species))
}

# ajout du taxID et des noms de chaque espèces 
matrix = cbind(TaxID, name, matrix)

head = c("Echantillon", "Program", "Species", "Percentage", "Rpkm", "Nb mapped reads")
echantillon = c()
species = c()
percentage = c()
rpkm_ech = c()
mapped_ech = c()

col = colnames(matrix)
       
for(i in 5:length(col)){
  sample = data.frame(matrix$name, matrix$TaxID, matrix[i])
  colnames(sample) = c("reference", "TaxID", "count")
  rpkm = (sample$count * (10**9))/(matrix$Size * as.numeric(sum(sample$count)))
  sample = cbind(sample, rpkm)
  sample = aggregate(cbind(rpkm, count) ~ TaxID + reference, sample, FUN=function(x) c(mn=mean(x), s=sum(x) ))
  sample = aggregate(cbind(rpkm[,1], count[,2]) ~ reference, sample, sum)
  perc = round(sample$V1/sum(sample$V1),3)
  sample = cbind(sample, perc)
  sample = sample[sample$perc != 0,]
 
#  sample_name = unlist(strsplit(col[i],"[.]|_"))[1]
  echantillon = c(echantillon,rep(col[i], dim(sample)[1]))
  species = c(species,as.character(sample$reference))
  percentage = c(percentage, sample$perc)
  rpkm_ech = c(rpkm_ech, sample$V1)
  mapped_ech = c(mapped_ech, sample$V2)
}

mod = rep(opt$label,length(echantillon))
df = data.frame(echantillon, mod, species, percentage, rpkm_ech, mapped_ech)
write.table(df, file=opt$out, row.names=FALSE, quote=FALSE, sep='\t', col.names = head)
