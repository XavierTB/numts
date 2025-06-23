##Program to detect numts in the COI Leray fragment
##It checks sequences against all genetic codes (or all metazoan codes if it is a metazoan) to detect stop codons
##For metazoans and seqs of 313 bp it also checks the five aminoacids conserved according to Pentinsaari et al 2016
##We assume that the codon starts in position 2 in the Leray fragment
##Most sequences should be 313 bp. Other lengths are acceptable if codon coherence is kept (i.e., length is 313+-3*n)

##You need to have these libreries installed

library(Biostrings)
library(ape)
library(ade4)



#Read a .csv with the ESVs. It should have a first column "id" with the ESV identification, a column "MOTU" with the id of the MOTU to which the corresponding ESV belong, and a final "seq" column with the sequence
#It will typically have columns with the number of reads found in each sample and whatever other columns with infos have been generated

ESVs <- read.csv("ASVs.csv",stringsAsFactors = F)
ESVs$MOTU<-factor(ESVs$MOTU,levels=unique(ESVs$MOTU))

#Read a .csv with the MOTUs. It should have a first column "id" with the the MOTU identification and a final "seq" column with the sequence
#It may also have columns with taxonomic identification. In particular, it needs to have a column "Kingdom" were metazoans are flagged "Metazoa"
#It must have a column "COUNT" with the total number of reads of the MOTU
#It will typically have other columns with the number of reads found in each sample
#It can have any other columns with descriptors

MOTUs <- read.csv("MOTUs.csv",stringsAsFactors = F)

##IMPORTANT MOTUs$id and ESVs$id should be coherent codes. For any MOTU, the MOTUs$id should be the ESV$id of the ESV representative of that MOTU. Do not use, for instance ESVxxxxxx for ESVs and MOTUxxxxxx for MOTUs


#Let's identify the columns with sample read information (here an instance, replace with the correct column assignments)
#This is needed only if you want the program to recalculate MOTU reads after deleting ESVs attributed to numts
col_samples_MOTUs<-c(16:282)#MOTUs has samples in columns 16:282
col_samples_ESVs<-c(4:270)#ESVs has samples in columns 4:270


#check
sum(MOTUs[,col_samples_MOTUs])==sum(ESVs[,col_samples_ESVs])#should be TRUE

###To make sure that the first column is named "id" and the last column is named "seq"
names(MOTUs)[1]<-"id"
names(ESVs)[1]<-"id"
names(MOTUs)[dim(MOTUs)[2]]<-"seq"
names(ESVs)[dim(ESVs)[2]]<-"seq"
###Put seqs in lower case:
ESVs$seq<-tolower(ESVs$seq)
MOTUs$seq<-tolower(MOTUs$seq)

#prepare files with results

numts_MOTUs<-as.data.frame(matrix(NA,length(levels(ESVs$MOTU)),17))

names(numts_MOTUs)<-c("MOTU_id","is_numt","gen_code","id_code","MOTUdepth_initial","ESVs_eliminated","MOTUdepth_final",
                    "ESVs_no_stop","ESVs_one_stop","ESVs_two_stops","ESVs_threeormore_stops","totalstops",
                    "is_metazoa_313","ESVs_wrong_aa","ESVs_wrong_aa_no_stops","wrong_aa_gen_code","wrong_aa_id_code")
numts_MOTUs$is_numt<-FALSE
numts_MOTUs$is_metazoa_313<-FALSE
numts_MOTUs[,5:12]<-0


numts_ESVs<-as.data.frame(matrix(NA,dim(ESVs)[1],10))
names(numts_ESVs)<-c("ESV_id","MOTU_ID","is_numt","gen_code","id_code","n_stops","is_metazoa_313","n_wrong_aa","wrong_aa_gen_code","wrong_aa_id_code")
numts_ESVs$is_numt<-FALSE
numts_ESVs$is_metazoa_313<-FALSE
numts_ESVs$n_stops<-0

##BEGIN NUMTS CORRECTION
for (i in 1:length(levels(ESVs$MOTU)))
{  
  datas<-ESVs[which(ESVs$MOTU==levels(ESVs$MOTU)[i]),]
  if (i/100-ceiling(i/100)==0) message(Sys.time()," processing MOTU ",i," of ",length(levels(ESVs$MOTU)))
  stops<-matrix(NA,dim(datas)[1],20)
  
  
  datas$seq<-as.character(datas$seq)  
  seq<-DNAStringSet(datas$seq)
  lenseq<-nchar(datas$seq)
  #We assume that codon starts in second position, as typical of the Leray fragment
  seq<-DNAStringSet(seq,start=2,end=lenseq)
  
  
  #Let's identify the representative ESV of the MOTU
  is_inlist<-FALSE
  representative<-NULL
  for (h in 1:dim(datas)[1])
    if (datas$seq[h]%in%MOTUs$seq) 
    {
      is_inlist<-TRUE 
      representative<-h
    }  
  
  if (is_inlist==FALSE) message("ERROR, no ESV of MOTU ",MOTUs$id[i],"has the representative seq of the MOTUs file")
  if (is_inlist)

#Let's check if it is a metazoan
  
is_metazoa<-FALSE
     

if (MOTUs$kingdom[which(MOTUs$seq==datas$seq[representative])]=="Metazoa")
          is_metazoa<-TRUE
  
  cod<-c(1:20)#for non-metazoans the 20 genetic codes are searched
  if (is_metazoa) cod<-c(1,2,4,5,7,11,12,15,18)#for metazoans only genetic codes of metazoans are searched
  for (qq in cod)
  {
    code<-getGeneticCode(as.character(GENETIC_CODE_TABLE$id[qq]))
    trans<-translate(seq,genetic.code=code)
    
    for (k in 1:dim(datas)[1])
    {
      aa<-strsplit(as.character(trans[k]),split="")
      aa<-unlist(aa)
      nstops<-length(which(aa=="*"))
      stops[k,qq]<-nstops
    }
  }
  
  
  bestcode<-which(colSums(stops)==min(colSums(stops),na.rm=T))[1] #we choose as best code the first one that gives the minimum mumber of stops
  bonscodes<-which(colSums(stops)==min(colSums(stops),na.rm=T))#to keep the other codes that produce the minimum of stops
  bestcodename<-GENETIC_CODE_TABLE$name[bestcode]
  bestcodeid<-GENETIC_CODE_TABLE$id[bestcode]
      
  #ONLY for metazoans and for seqs of 313 bp we check for bad aminoacids in the five positions conserved according to Pentinsaari et al., 2016
  #For assigning bad aa we check not only the first best code, but all codes that produced the minimum number of stop codons (bonscodes)
  #We apply the same code to all ESVs of the MOTU. We keep the code generating the mininum number of bad aa over all seqs of the MOTU. 
    
  aa_xung<-matrix(NA,dim(datas)[1],length(bonscodes))
  if (is_metazoa & nchar(datas$seq[representative])==313)
      for (qq in 1:length(bonscodes))
      {
        code<-getGeneticCode(as.character(GENETIC_CODE_TABLE$id[bonscodes[qq]]))
        trans<-translate(seq,genetic.code=code)
        for (k in 1:dim(datas)[1])
        {
          aa<-strsplit(as.character(trans[k]),split="")
          aa<-unlist(aa)
          bad_aa<-0
          if (aa[20]!="H") bad_aa<-bad_aa+1
          if (aa[23]!="G") bad_aa<-bad_aa+1
          if (aa[32]!="N") bad_aa<-bad_aa+1
          if (aa[81]!="D") bad_aa<-bad_aa+1
          if (aa[95]!="G") bad_aa<-bad_aa+1
          aa_xung[k,qq]<-bad_aa
        }
      }
      
      bestindex_aa<-which(colSums(aa_xung)==min(colSums(aa_xung)))[1]
      bestcode_aa<-bonscodes[which(colSums(aa_xung)==min(colSums(aa_xung)))[1]]
      bestcodename_aa<-GENETIC_CODE_TABLE$name[bestcode_aa]
      bestcodeid_aa<-GENETIC_CODE_TABLE$id[bestcode_aa]
  
  
  
  #let's fill the results tables
  
  
  numts_MOTUs$MOTU_id[i]<-datas$id[representative]
  numts_MOTUs$MOTUdepth_initial[i]<-dim(datas)[1]
  numts_MOTUs$gen_code[i]<-bestcodename
  numts_MOTUs$id_code[i]<-bestcodeid
  numts_MOTUs$ESVs_no_stop[i]<-sum(stops[,bestcode]==0)
  numts_MOTUs$ESVs_one_stop[i]<-sum(stops[,bestcode]==1)
  numts_MOTUs$ESVs_two_stops[i]<-sum(stops[,bestcode]==2)
  numts_MOTUs$ESVs_threeormore_stops[i]<-sum(stops[,bestcode]>2)
  numts_MOTUs$totalstops[i]<-sum(stops[,bestcode])
  if (is_metazoa & nchar(datas$seq[representative])==313) {
    numts_MOTUs$is_metazoa_313[i]<-TRUE
    numts_MOTUs$ESVs_wrong_aa[i]<-length(which(aa_xung[,bestindex_aa]>0))
    numts_MOTUs$ESVs_wrong_aa_no_stops[i]<-length(which(stops[,bestcode]==0 & aa_xung[,bestindex_aa]>0))
    numts_MOTUs$wrong_aa_gen_code[i]<-bestcodename_aa
    numts_MOTUs$wrong_aa_id_code[i]<-bestcodeid_aa
  } 
  
  numts_MOTUs$MOTUdepth_final[i]<-sum(numts_MOTUs$ESVs_no_stop[i],-numts_MOTUs$ESVs_wrong_aa_no_stops[i],na.rm=T)
  numts_MOTUs$ESVs_eliminated[i]<-numts_MOTUs$MOTUdepth_initial[i]-numts_MOTUs$MOTUdepth_final[i]
  if (numts_MOTUs$ESVs_eliminated[i]==numts_MOTUs$MOTUdepth_initial[i]) numts_MOTUs$is_numt[i]<-TRUE
positions<-which(ESVs$MOTU==levels(ESVs$MOTU)[i])
flag<-vector("logical",length(positions))
for (kk in 1:length(positions))
{
pos<-positions[kk]
numts_ESVs$ESV_id[pos]<-ESVs$id[pos]
numts_ESVs$MOTU_id[pos]<-datas$id[representative]
numts_ESVs$gen_code[pos]<-bestcodename
numts_ESVs$id_code[pos]<-bestcodeid
numts_ESVs$n_stops[pos]<-stops[kk,bestcode]
numts_ESVs$is_metazoa_313[pos]<-(is_metazoa & nchar(datas$seq[representative])==313)
if (is_metazoa & nchar(datas$seq[representative])==313) {
  numts_ESVs$n_wrong_aa[pos]<-aa_xung[kk,bestindex_aa]
  numts_ESVs$wrong_aa_gen_code[pos]<-bestcodename_aa
  numts_ESVs$wrong_aa_id_code[pos]<-bestcodeid_aa
}
flag[kk]<-sum(numts_ESVs$n_stops[pos],numts_ESVs$n_wrong_aa[pos],na.rm=T)>0
if (flag[kk]) numts_ESVs$is_numt[pos]<-TRUE
} 
}
#END OF NUMTS CORRECTION


#delete MOTUs with errors (likely numts). These are the MOTUs for which all ESVs have errors    
flag<-vector("logical",dim(numts_MOTUs)[1])
for (i in 1:dim(numts_MOTUs)[1]) 
{
  MOT<-numts_MOTUs[which(numts_MOTUs$MOTU_id==MOTUs$id[i]),]
  if (MOT$is_numt) flag[i]<-TRUE
}
MOTUss<-MOTUs
if (sum(flag)>0) MOTUss<-MOTUss[-which(flag),]

#check:
dim(MOTUss)[1]==sum(!numts_MOTUs$is_numt)#should be TRUE

#delete ESVs with errors   

flag<-vector("logical",dim(numts_ESVs)[1])
for (i in 1:dim(numts_ESVs)[1]) 
{
  ESV<-numts_ESVs[which(numts_ESVs$ESV_id==ESVs$id[i]),]
  if (ESV$is_numt) flag[i]<-TRUE
}
ESVss<-ESVs
if (sum(flag)>0) ESVss<-ESVss[-which(flag),]
ESVss$MOTU<-factor(ESVss$MOTU,levels=unique(ESVss$MOTU))

#check
dim(ESVss)[1]==sum(!numts_ESVs$is_numt)#should be TRUE


#From the ESV file we reconstruct the MOTUs file by recounting abundances 

for (i in 1:length(levels(ESVss$MOTU)))
{
    MOTUss[which(MOTUss$id==levels(ESVss$MOTU)[i]),col_samples_MOTUs]<-
          colSums(ESVss[which(ESVss$MOTU==levels(ESVss$MOTU)[i]),col_samples_ESVs])
    MOTUss$CLUST_WEIGHT[i]<-length(which(ESVss$MOTU==levels(ESVss$MOTU)[i]))
    }

#redo the column COUNT
  MOTUss$COUNT<-rowSums(MOTUss[,col_samples_MOTUs])
  
##If the representative sequence of a MOTU is eliminated, we eliminate the MOTU even if someof its ESVs are deemed as valid
##never happened so far, but it is possible

#representative sequences eliminated will exist in MOTUss$seq but not in ESVss$seq)
dolentes<-which(!MOTUss$seq%in%ESVss$seq) 
if (length(dolentes)>0) 
{
MOTUss<-MOTUss[,-dolentes]
ESVss<-ESVss[,-which(ESVss$seq%in%MOTUss$seq[dolentes])]
}


#check
  dim(MOTUss)[1]==length(levels(ESVss$MOTU))#should be TRUE
  which(!MOTUss$id%in%levels(ESVss$MOTU))#should be 0
  
#More checks, all these must produce the same result
  sum(MOTUss[,col_samples_MOTUs])
  sum(MOTUss$COUNT)
  sum(ESVss[,col_samples_ESVs])
  sum(ESVss$COUNT)
 
#More checks, all these must produce the same result
  dim(ESVss)[1]
  sum(MOTUss$CLUST_WEIGHT)
  length(which(!numts_ESVs$is_numt))


#write results

write.csv(ESVss,"ASVs_nonumts.csv",row.names = F)    
write.csv(MOTUss,"MOTUs_nonumts.csv",row.names = F)   
write.csv(numts_MOTUs,"numts_MOTUs.csv",row.names = F)
write.csv(numts_ESVs,"numts_ASVs.csv",row.names = F)

print(paste("we deleted ",sum(numts_MOTUs$is_numt)," motus"))
print(paste("we deleted ",sum(numts_ESVs$is_numt)," ESVs"))

print(paste("ESVs deleted because of stops ",sum(numts_ESVs$n_stops>0,na.rm=T)))
print(paste("total stops ",sum(numts_ESVs$n_stops,na.rm=T)))
print(paste("ESVs deleted because of wrong_aa without stops ",sum(numts_ESVs$n_stops==0&numts_ESVs$n_wrong_aa>0,na.rm=T)))
print(paste("ESVs with wrong_aa and stops ",sum(numts_ESVs$n_stops>0&numts_ESVs$n_wrong_aa>0,na.rm=T)))
print(paste("total ESVs with wrong_aa ",sum(numts_ESVs$n_wrong_aa>0,na.rm=T)))
print(paste("total wrong_aa ",sum(numts_ESVs$n_wrong_aa,na.rm=T)))
print(paste("MOTUs deleted because representative sequence was a numt ",length(dolentes)))

message("end")

