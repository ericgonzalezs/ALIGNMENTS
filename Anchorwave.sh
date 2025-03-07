#REF PREP

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --cpus-per-task=10
#SBATCH --mem=90G
module load StdEnv/2020 minimap2/2.24 gcc/11.3.0

export PATH=$PATH:/home/egonza02/scratch/SOFTWARE/ANCHORWAVE/anchorwave

anchorwave gff2seq -r ANN1372_HAP1_hap1.reviewed.chr_assembled.fasta -i ANN1372H1_mincov90_minID90.PANNEW76k.gmap.gff3 -o refH1.cds.fa

anchorwave gff2seq -r ANN1372_HAP2_hap2.reviewed.chr_assembled.fasta -i ANN1372H2_mincov90_minID90.PANNEW76k.gmap.gff3 -o refH2.cds.fa


minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 ANN1372_HAP1_hap1.reviewed.chr_assembled.fasta refH1.cds.fa > refH1.sam

minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 ANN1372_HAP2_hap2.reviewed.chr_assembled.fasta refH2.cds.fa > refH2.sam

##################################
#Alignment

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=2-0
#SBATCH --cpus-per-task=10
#SBATCH --mem=497G
module load StdEnv/2020 minimap2/2.24 gcc/11.3.0
export PATH=$PATH:/home/egonza02/scratch/SOFTWARE/ANCHORWAVE/anchorwave
export TMPDIR=/home/egonza02/scratch/ALIGNMENTS/GMAP/ALIGNMENTS

i=$(cat SEPHAP_GENOMES.txt | grep -v "ANN1372_V3_corrected_hap1.reviewed.chr_assembled.fasta" | head -n 1 | tail -n 1)
name=$(echo $i | cut -d "." -f 1-3)

minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 $i refH1.cds.fa >  "$name"".sam"


anchorwave proali -i ANN1372H1_V3_corrected_mincov90_minID90.PANNEW76k.gmap.gff3 -r ANN1372_V3_corrected_hap1.reviewed.chr_assembled.fasta -a "$name"".sam" -as refH1.cds.fa -ar refH1.sam -s $i -n "$name""_vs_ANN1372_HAP1.anchors" -o "$name""_vs_ANN1372_HAP1_anchorwave.maf" -t 9 -R 1 -Q 1 -f "$name""_vs_ANN1372_HAP1_anchorwave.f.maf"  >  "$name""_vs_ANN1372_HAP1_anchorwave.log"

#######################################
#convert maf to table
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --mem=50G

module load StdEnv/2020 python/2.7.18
i=$(ls *anchorwave.maf | head -n 1 | tail -n 1)

name=$(echo $i | cut -d "." -f 1)

python2 /home/egonza02/scratch/SOFTWARE/ANCHORWAVE/anchorwave/scripts/maf-convert tab $i  > "$name"".tab"
#you have to download maf-convert from here https://gitlab.com/mcfrith/last/-/blob/main/bin/maf-convert

#################################################
#convert table to the format to plot
for i in $(ls *anchorwave.tab)
do
       name=$(echo $i | cut -d "." -f 1 )

       grep -v "#" $i | awk -F "\t" -v OFS="\t" '{print $3, $3 + $4, $8, $8 + $9, 30, $7, $2, $10}' > "$name""_tab_fp.txt"

done


##################################################
#PLOT RUN_PLOTS_Per_Chr_allchr.R

library(dplyr)
library(magrittr)
library(GenomicRanges)
library(knitr)
library(ggplot2)
library(tidyr)
library(tidyverse)

#functions
readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

filterMum <- function(df, minl=1000, flanks=1e4){
    coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
        summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
        ungroup %>% arrange(desc(rs)) %>%
        mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
    merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
        mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}

diagMum <- function(df){
    ## Find best qid order
    rid.o = df %>% group_by(qid, rid) %>% summarize(base=sum(abs(qe-qs)),
                                                    rs=weighted.mean(rs, abs(qe-qs))) %>%
        ungroup %>% arrange(desc(base)) %>% group_by(qid) %>% do(head(., 1)) %>%
        ungroup %>% arrange(desc(rid), desc(rs)) %>%
        mutate(qid=factor(qid, levels=unique(qid)))
    ## Find best qid strand
    major.strand = df %>% group_by(qid) %>%
        summarize(major.strand=ifelse(sum(sign(qe-qs)*abs(qe-qs))>0, '+', '-'),
                  maxQ=max(c(qe, qs)))
    merge(df, major.strand) %>% mutate(qs=ifelse(major.strand=='-', maxQ-qs, qs),
                                       qe=ifelse(major.strand=='-', maxQ-qe, qe),
                                       qid=factor(qid, levels=levels(rid.o$qid)))
}


        args <- commandArgs(trailingOnly = TRUE)

                f <-  args[1] #delta
               chr <- args[2]
                #bed <- args[3] #trans
        #        bed2 <-  args[3] #inv
        #       bed3 <- args[4] # centromers
        #       bed4 <- args[5] #chr centromers
                #st <- args[3]
                #end <- args[4]

         #      centromers <-  read.table(bed3,  header=F)
          #     centromers <- centromers[order(centromers[,1], centromers[,2]), ]
          #     centromers$start <- pmin(centromers$V2, centromers$V3)
          #     centromers$end <- pmax(centromers$V2, centromers$V3)
          #     centromers <- subset(centromers, centromers$V1== bed4)
          #     centromers_s <- centromers$start
               #TRANS <- read.table(bed,  header=F)
               #TRANS_L <- TRANS$V7
               #TRANS_R <- TRANS$V8

         #      INV <- read.table(bed2,  header=F)
         #      INV <- subset(INV, INV$V7 == chr)
         #      INV$lenght <- INV$V2 - INV$V1
         #      INV_long <- subset(INV, INV$lenght > 100000)

         #      INV_L <- INV_long$V1
         #      INV_R <- INV_long$V2
               #INV_L <- INV$V7
               #INV_R <- INV$V8

                 name <- gsub(".txt", "", f)

                # name2 <-  paste(name, chr, "Ha412", ".jpg" , sep="_")
                 #mumgp = readDelta(f)
                 mumgp <- read.table(f, header=F)
                 colnames(mumgp) <- c("rs", "re", "qs", "qe", "error", "qid", "rid", "strand")
                 mumgp[c("qs", "qe")] <- t(mapply(\(a, b, c, d, e, f, g, h){
                 if(h == "-") c(d,c) else c(c,d) }, mumgp$rs, mumgp$re, mumgp$qs, mumgp$qe, mumgp$error, mumgp$qid, mumgp$rid, mumgp$strand))
                 #mumgp$qid <- paste(mumgp$qid, "Arg", sep="_")

                 #mumgp <- read.table(f, header=F)
                 #colnames(mumgp) <- c("rs", "re", "qs", "qe", "error", "qid", "rid", "strand")
                 #mumgp[c("qs", "qe")] <- t(mapply(\(a, b, c, d, e, f, g, h){
                  #  if(h == "-") c(d,c) else c(c,d) }, mumgp$rs, mumgp$re, mumgp$qs, mumgp$qe, mumgp$error, mumgp$qid, mumgp$rid, mumgp$strand))

                 mumgp <- subset(mumgp, mumgp$rid == chr) #& mumgp$rs > st & mumgp$re < end ) FILTER PER CHR

                 mumgp$H1_H1_name <- paste(sub("^(H\\d+).*", "\\1", mumgp$rid), round(mumgp$rs/1000000, digits = 2), round(mumgp$re/1000000, digit
s = 2), sub("^(H\\d+).*", "\\1", mumgp$qid), round(mumgp$qs/1000000, digits = 2), round(mumgp$qe/1000000, digits = 2), sep = "-")

                 mumgp.filt = filterMum(mumgp, minl=1e4)
                 mumgp.filt.diag = diagMum(mumgp.filt)

 #  mumgp.filt.diag$H1_H1_name <- paste(sub("^(H\\d+).*", "\\1", mumgp.filt.diag$rid), round(mumgp.filt.diag$rs/1000000, digits = 1), round(mumgp.
filt.diag$re/1000000, digits = 1), sub("^(H\\d+).*", "\\1", mumgp.filt.diag$qid), round(mumgp.filt.diag$qs/1000000, digits = 1), round(mumgp.filt.
diag$qe/1000000, digits = 1), sep = "-")


        #add vertical lines to the transposition region
            #  data_vline <- data.frame(group = c("Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr6", "Chr6", "Chr6", "Chr8", "Chr9", "Chr11", "Chr12"
, "Chr14", "Chr15", "Chr16", "Chr17"), vline
#=c(NA, NA, NA, NA, 125000000, 130000000, 130000000, 155000000, NA, NA, NA, NA, NA, NA, NA, NA))

              P1 <- ggplot(mumgp.filt.diag, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
  geom_segment(show.legend=FALSE, size=3) + geom_point(alpha=0.09) + theme_bw() +
   geom_text(aes(x = re, y = qe, label = H1_H1_name), size = 3) +
  facet_grid(qid~rid, scales='free', space='free', switch='both') +
  guides(colour=guide_legend(override.aes=list(alpha=1))) +
  theme(strip.text.y=element_text(angle=180, size=10),
        strip.text.x=element_text(size=10),
        strip.background=element_blank(),
        legend.position=c(1,-.03), legend.justification=c(1,1),
        legend.direction='horizontal',
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.spacing=unit(0, 'cm')) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1') #+
        # geom_vline(xintercept=INV_R, color = "green", size=1) +
        # geom_vline(xintercept=INV_L, color = "green", size=1) +
#        geom_vline(xintercept=centromers_s, color = "purple", size=1, alpha = .05) +
        # annotate("rect", ymin = 0, ymax =  max( mumgp.filt.diag$qe), xmin = INV_L, xmax = INV_R, alpha = 0.2, fill = "green")
        # geom_vline(xintercept=end, color = "purple", size=1, alpha = .05) +
        #geom_vline(data = mumgp.filt.diag %>% filter(rid == "Ha412HOChr03"), aes(xintercept = 28000000), linetype="dotted", size=2) +
        #geom_vline(data = mumgp.filt.diag %>% filter(rid == "Ha412HOChr03"), aes(xintercept = 42000000), linetype="dotted", size=2) +
        #geom_vline(data = mumgp.filt.diag %>% filter(rid == "Ha412HOChr03"), aes(xintercept = 89000000), linetype="dotted", size=2)
       # geom_hline(yintercept=TRANS_L, color = "green", size=1,  alpha = .05) +
       # geom_hline(yintercept=TRANS_R, color = "green", size=1,  alpha = .05) +
       # geom_hline(yintercept=INV_R, color = "red", size=1, alpha = .2) +
       # geom_hline(yintercept=INV_L, color = "red", size=1, alpha = .2) +
       # annotate("rect", xmin = 0, xmax =  max( mumgp.filt.diag$re), ymin = INV_L, ymax = INV_R, alpha = 0.2, fill = "red")

jpeg( paste(name, chr, ".jpg", sep ="_"),  width=2000, height= 2000)
print(P1)
dev.off()


###############################################################
#How to run it

module load StdEnv/2020  r/4.1.2
for j in $( ls *tab_fp.txt)
do

   for i in $(cut -f 7 $j |sort |  uniq)
   do


   Rscript RUN_PLOTS_Per_Chr_allchr.R $j $i

   done

done
