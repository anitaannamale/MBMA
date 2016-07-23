
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type = "character", default=NULL, 
              help="count matrice file name", metavar="character"),
  make_option(c("-l","--label"), type = "character", default = "MBMA",
              help = "label for the analysis [default= %default]", metavar = "character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

# Lecture des fichiers
matrix = read.table(opt$file, header=T ,sep='\t', stringsAsFactors = F)

# Récupération du taxID des souches
geneID = NULL
for(id in matrix$ref_genome) {
  geneID = c(geneID, unlist(strsplit(id, "_"))[1])
}

# ajout du taxID et des noms de chaque espèces 
matrix = cbind(geneID, matrix)

head = c("Echantillon", "Program", "Resistance genes", "Percentage", "Rpkm", "Nb mapped reads")
echantillon = c()
resistance = c()
percentage = c()
rpkm_ech = c()
mapped_ech = c()

col = colnames(matrix)
       
for(i in 4:length(col)){
  sample = data.frame(matrix$geneID, matrix[i])
  colnames(sample) = c("geneID", "count")
  rpkm = (sample$count * (10**9))/(matrix$Size * as.numeric(sum(sample$count)))
  sample = cbind(sample, rpkm)
  sample = aggregate(cbind(rpkm, count) ~ geneID, sample, FUN=sum)
  perc = round(sample$rpkm/sum(sample$rpkm),3)
  sample = cbind(sample, perc)
  sample = sample[sample$perc != 0,]
 
#  sample_name = unlist(strsplit(col[i],"[.]|_"))[1]
  echantillon = c(echantillon,rep(col[i], dim(sample)[1]))
  resistance = c(resistance,as.character(sample$geneID))
  percentage = c(percentage, sample$perc)
  rpkm_ech = c(rpkm_ech, sample$rpkm)
  mapped_ech = c(mapped_ech, sample$count)
}

mod = rep(opt$label,length(echantillon))
df = data.frame(echantillon, mod, resistance, percentage, rpkm_ech, mapped_ech)
write.table(df, file=opt$out, row.names=FALSE, quote=FALSE, sep='\t', col.names = head)
