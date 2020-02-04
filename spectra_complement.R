#20200203
#Introduce variables
#20190725
#We optimized the in-silico library:
#1. lib containing y-ions only
#2. all missed-cleaved peptides were included, heavy and light.
args=commandArgs(T)
library(stringr)
options(stringsAsFactors = F)

#variables
spectrallibrary <- c('C:/Users/lenovo/Desktop/tlib.tsv') #args[1]
outputlibrary <- args[2]

#read spectal library file
L0 <- read.table(spectrallibrary, header = T, sep="\t", fill =T)
L0 <- L0[!grepl('DECOY', L0$transition_name), ]


#add heavy m/z to y ions
L0_f00 <- L0
L0_f00up <- L0_f00
FragmentIonType <- vector()
FragmentSeriesNumber <- vector()
for (i in 1:nrow(L0_f00)){
  ProductCharge <- unlist(strsplit(L0_f00$transition_name[i], split = '_', fixed = F))[3]
  IonInfo <- unlist(strsplit(L0_f00$transition_name[i], split = '_', fixed = F))[2]
  FragmentIonType <- unlist(strsplit(IonInfo, split = '', fixed = F))[1]
  FragmentSeriesNumber <- unlist(strsplit(IonInfo, split = '', fixed = F))[2]
  
  sep_seq <- str_sub(L0_f00$PeptideSequence[i], 1, as.numeric(FragmentSeriesNumber))
  sep_seq_rev <- str_sub(L0_f00$PeptideSequence[i], -as.numeric(FragmentSeriesNumber), -1)                
  kax <- str_count(sep_seq, pattern = 'K') + str_count(sep_seq, pattern = 'R')
  kax_rev <- str_count(sep_seq_rev, pattern = 'K') + str_count(sep_seq_rev, pattern = 'R')

  L0_f00up$PrecursorMz[i] <- L0_f00$PrecursorMz[i] + 6.02012898*(str_count(L0_f00$PeptideSequence[i], pattern = 'K') + str_count(L0_f00$PeptideSequence[i], pattern = 'R'))/L0_f00$PrecursorCharge[i]
  
  if(FragmentIonType == 'b'){
    L0_f00up$ProductMz[i] <- L0_f00$ProductMz[i] + 6.02012898*(kax)/as.numeric(ProductCharge)
  }else{
    L0_f00up$ProductMz[i] <- L0_f00$ProductMz[i] + 6.02012898*(kax_rev)/as.numeric(ProductCharge)
  }
  if(i%%100 == 0) print(i)
}

#add unimod labels
L0_f00up$transition_name <- gsub('K', 'K(UniMod:188)', L0_f00up$transition_name, fixed = T)
L0_f00up$transition_name <- gsub('R', 'R(UniMod:188)', L0_f00up$transition_name, fixed = T)
L0_f00up$transition_group_id <- gsub('K', 'K(UniMod:188)', L0_f00up$transition_group_id, fixed = T)
L0_f00up$transition_group_id <- gsub('R', 'R(UniMod:188)', L0_f00up$transition_group_id, fixed = T)
L0_f00up$FullUniModPeptideName <- gsub('K', 'K(UniMod:188)', L0_f00up$FullUniModPeptideName, fixed = T)
L0_f00up$FullUniModPeptideName <- gsub('R', 'R(UniMod:188)', L0_f00up$FullUniModPeptideName, fixed = T)
L0_f00up$PeptideGroupLabel <- gsub('K', 'K(UniMod:188)', L0_f00up$PeptideGroupLabel, fixed = T)
L0_f00up$PeptideGroupLabel <- gsub('R', 'R(UniMod:188)', L0_f00up$PeptideGroupLabel, fixed = T)

#combine to generate a heavy library
H1 <- L0_f00up
H1$LabelType = 'H'
H1_check <- H1[!grepl('(UniMod:188)', H1$transition_name), ]
H1_f1 <- H1[grepl('(UniMod:188)', H1$transition_name), ]
H1_f2 <- H1_f1[!grepl('iRT', H1_f1$ProteinName, fixed = T), ]
iRT <- L0_f00[grepl('iRT', L0_f00$ProteinName, fixed = T), ]
H1_f3 <- rbind(H1_f2, iRT)
H1_f3$LabelType = 'H'
#write.table(H1_f3, file = 'PX20191011_K562_MQLib_pseudoH.tsv', sep = '\t', row.names = F, col.names = T, quote = F)

#combine to generate a consensus library
#delete iRTs in heavy lib
H1_noiRT <- H1_f3[!grepl('iRT', H1_f3$ProteinName), ]
H1_noiRT <- na.omit(H1_noiRT)
L0$LabelType = 'L'
LH1 <- data.frame(rbind(as.matrix(H1_noiRT), as.matrix(L0)))
write.table(LH1, file = outputlibrary, sep = '\t', row.names = F, col.names = T, quote = F)
print('spectra_complement done.', quote = F)