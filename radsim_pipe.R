# Get fragment of RAD seq

#libraries
library(ape)
library(SimRAD)
library(stringr)
library(seqinr)

#loading sequences from a folder
filelocs <- list.files("/home/fgrewe/lichen_Rhizo/rhizoplaca_consensus_mappings", pattern=".fasta", full.names=TRUE)
iters <- length(filelocs)
filenames <- str_sub((list.files("/home/fgrewe/lichen_Rhizo/rhizoplaca_consensus_mappings", pattern=".fasta", full.names=FALSE)), 1, -7)

#adding enzyme sequence ApeK1
cs_5p1 <- "G"
cs_3p1 <- "CAGC"
cs_5p2 <- "G"
cs_3p2 <- "CTGC"

#going through the loop
for (i in 1:iters){
    
    print(filenames[i])
        
    #loading sequence
    refseq <- ref.DNAseq(filelocs[i], subselect.contigs = FALSE)
    
    #digest and size selection
    simseq.dig <- insilico.digest(refseq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose = TRUE)
    size.superwide <- size.select(simseq.dig, 200, 500, graph = FALSE, verbose = TRUE)
    
	#list of all fragements smaller than 243
	size.smaller <- size.superwide[nchar(size.superwide) < 243]
	size.larger <- size.superwide[nchar(size.superwide) >= 243]
	
	#get full frgament when overlap less then 30% (replicate overlapping clusters from pyRAD: " I found that reads which ovelap by less than about 30% are difficult to detect, and the minimum amount of overlap allowed is not a parameter you can change in the params file, but it can be edited in the pyrad code.")
	sub.super.merged <- size.smaller
	
	#get subset
	sub.super.start <- str_sub(size.larger, 1, 143)
	sub.super.stop <- str_sub(size.larger, -143, -1)
	
	
	#write first and last 50bp in fasta
	write.fasta(as.list(sub.super.start), paste(filenames[i], c((1:length(sub.super.start)+100000)), "r1", sep="_"), file.out=paste(filenames[i], "R1.fas", sep="_"))
	write.fasta(as.list(sub.super.stop), paste(filenames[i], c((1:length(sub.super.stop)+100000)), "r2", sep="_"), file.out=paste(filenames[i], "R2.fas", sep="_"))
	write.fasta(as.list(sub.super.merged), paste(filenames[i], c((1:length(sub.super.merged)+100000)), "m", sep="_"), file.out=paste(filenames[i], "M.fas", sep="_"))
	#write.fasta(as.list(size.superwide), paste(filenames[i], c(1:length(sub.super.stop)), "frag", sep="_"), file.out=paste(filenames[i], "ApeK1_fragments.fasta", sep="_"))  

	system(paste("cat ",filenames[i],"_R1.fas >> ",filenames[i],".consens", sep=""))
	system(paste("cat ",filenames[i],"_R2.fas >> ",filenames[i],".consens", sep=""))
	system(paste("cat ",filenames[i],"_M.fas >> ",filenames[i],".consens", sep=""))
	system(paste("gzip ",filenames[i],".consens", sep=""))
	#further processing not required when pluggin after step 5
	#system(paste("perl fasta_to_fastq.pl ",filenames[i],".fasta > ",filenames[i],".fastq", sep=""))
	#system(paste("gzip ",filenames[i],".fastq", sep=""))
}

# ##################
# # In command line:
# ##################


# #system(paste("cat ",filenames[3],"_R1.fas >> ",filenames[3],"_both.fas", sep=""))
# system("gzip *.consens")
# system("mkdir ./clust.84")
# system("cp *consens.gz ./clust.84")
# system("mkdir ./stats")
# system("pyRAD -p params.txt -s 67")
# #runs in error when writing .geno file

# #####################################	
# #####################################
# ## for a single sequence start here:
# #####################################
	
# library(SimRAD)
# library(stringr)
# library(seqinr)

# #loading sequence
# refseq <- ref.DNAseq("./sequences/mela_REF.fasta", subselect.contigs = FALSE)

# #adding enzyme sequence ApeK1
# cs_5p1 <- "G"
# cs_3p1 <- "CAGC"
# cs_5p2 <- "G"
# cs_3p2 <- "CTGC"

# #addting enzyme pstI for testing purpose
# #cs_5p1 <- "CTGCA"
# #cs_3p1 <- "G"
# #simseq.dig <- insilico.digest(refseq, cs_5p1, cs_3p1, verbose = TRUE)

# #digest and size selection
# simseq.dig <- insilico.digest(refseq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose = TRUE)
# size.superwide <- size.select(simseq.dig, 200, 500, graph = TRUE, verbose = TRUE)

# #get subset
# sub.super.start <- str_sub(size.superwide, 1, 143)
# sub.super.stop <- str_sub(size.superwide, -250, -1)

# #write first and last 50bp in fasta
# write.fasta(as.list(sub.super.start), paste("refseq_first", c(1:length(sub.super.start)), sep="_"), file.out="ref250first.fasta")
# write.fasta(as.list(sub.super.stop), paste("refseq_last", c(1:length(sub.super.stop)), sep="_"), file.out="ref250last.fasta")
# write.fasta(as.list(size.superwide), paste("refseq_frag", c(1:length(sub.super.stop)), sep="_"), file.out="PstIfragments.fasta")

# ##########################################################
# ##########################################################
# # parameter file for .84 cluster and min 4 taxa per locus
# ##########################################################


# ==** parameter inputs for pyRAD version 3.0.64  **======================== affected step ==
# ./                        ## 1. Working directory                                 (all)
# ../*.fastq              ## 2. Loc. of non-demultiplexed files (if not line 18)  (s1)
# ../*.barcodes              ## 3. Loc. of barcode file (if not line 18)             (s1)
# vsearch                   ## 4. command (or path) to call vsearch (or usearch)    (s3,s6)
# muscle                    ## 5. command (or path) to call muscle                  (s3,s7)
# GWGC                     ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG)     (s1,s2)
# 12                     ## 7. N processors (parallel)                           (all)
# 1                         ## 8. Mindepth: min coverage for a cluster              (s4,s5)
# 4                         ## 9. NQual: max # sites with qual < 20 (or see line 20)(s2)
# .84                    ## 10. Wclust: clustering threshold as a decimal        (s3,s6)
# rad                       ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(others:see docs)(all)
# 4                         ## 12. MinCov: min samples in a final locus             (s7)
# 3                         ## 13. MaxSH: max inds with shared hetero site          (s7)
# c84d6m4p3                 ## 14. Prefix name for final output (no spaces)         (s7)
# ==== optional params below this line ===================================  affected step ==
                       # ## 15.opt.: select subset (prefix* only selector)            (s2-s7)
                       # ## 16.opt.: add-on (outgroup) taxa (list or prefix*)         (s6,s7)
                       # ## 17.opt.: exclude taxa (list or prefix*)                   (s7)
                       # ## 18.opt.: loc. of de-multiplexed data                      (s2)
                       # ## 19.opt.: maxM: N mismatches in barcodes (def= 1)          (s1)
                       # ## 20.opt.: phred Qscore offset (def= 33)                    (s2)
                       # ## 21.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict   (s2)
                       # ## 22.opt.: a priori E,H (def= 0.001,0.01, if not estimated) (s5)
                       # ## 23.opt.: maxN: max Ns in a cons seq (def=5)               (s5)
                       # ## 24.opt.: maxH: max heterozyg. sites in cons seq (def=5)   (s5)
                       # ## 25.opt.: ploidy: max alleles in cons seq (def=2;see docs) (s4,s5)
                       # ## 26.opt.: maxSNPs: (def=100). Paired (def=100,100)         (s7)
                       # ## 27.opt.: maxIndels: within-clust,across-clust (def. 3,99) (s3,s7)
                       # ## 28.opt.: random number seed (def. 112233)              (s3,s6,s7)
                       # ## 29.opt.: trim overhang left,right on final loci, def(0,0) (s7)
# n,p,s,v              ## 30.opt.: output formats: p,n,a,s,v,u,t,m,k,g,* (see docs) (s7)
                       # ## 31.opt.: maj. base call at depth>x<mindepth (def.x=mindepth) (s5)
                       # ## 32.opt.: keep trimmed reads (def=0). Enter min length.    (s2)
                       # ## 33.opt.: max stack size (int), def= max(500,mean+2*SD)    (s3)
                       # ## 34.opt.: minDerep: exclude dereps with <= N copies, def=1 (s3)
                       # ## 35.opt.: use hierarchical clustering (def.=0, 1=yes)      (s6)
                       # ## 36.opt.: repeat masking (def.=1='dust' method, 0=no)      (s3,s6)
                       # ## 37.opt.: vsearch max threads per job (def.=6; see docs)   (s3,s6)
# ==== optional: list group/clade assignments below this line (see docs) ==================

	

