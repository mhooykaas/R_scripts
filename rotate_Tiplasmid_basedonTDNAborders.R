# to locate  T-DNA borders in Ti/Ri plasmid and rotate the molecule to make it start at the left T-DNA border of the first (or only) T-DNA
# also locate overdrive sequence (should be close to T-DNA border)
# for interactive use in Rstudio
library(seqinr)
library(Biostrings)

setwd("D:/working/directory")
# read plasmid sequence:
ri <- readDNAStringSet("D:/location/file/to/be/rotated/contig4.fasta",format="fasta")

# define sequences to be searched
left_TDNAborder <- DNAString("NGGCAGGATATATNNNNNTGTAAAN")
right_TDNAborder <- DNAString("TGNCAGGATATATNNNNNNGTNNNN")
consensus_TDNAborder <- DNAString("NGNCAGGATNTATNNNNNNGTNNNN") #TNTAT instead of TATAT
consensus_TDNAborder_rc <- reverseComplement(consensus_TDNAborder)
overdrive <- DNAString("TNNNNNNNTNNNTNTGTTT")
overdrive_rc <- reverseComplement(overdrive)

# print locations of these patterns to screen
vmatchPattern(left_TDNAborder,ri,fixed=FALSE)
vmatchPattern(right_TDNAborder,ri,fixed=FALSE)
vmatchPattern(consensus_TDNAborder,ri,fixed=FALSE)
vmatchPattern(consensus_TDNAborder_rc,ri,fixed=FALSE)
vmatchPattern(overdrive,ri,fixed=FALSE)
vmatchPattern(overdrive_rc,ri,fixed=FALSE)

#rearrange Ri plasmid
width(ri)
end <- subseq(ri,start=74624,end=252168)  # end segment = from where T-DNA starts to the end of the molecule
beginning <- subseq(ri,start=1,end=74623) # start segment = from position 1 to T-DNA start site minus 1

rearranged <- paste0(end,beginning)
width(rearranged) # should be equal to width(ri)

write.fasta(rearranged,"contig4_rotated","contig4_Ti_rotated.fasta")

#check resulting plasmid: first border should be found at 1
ti <- readDNAStringSet("D:/working/directory/contig4_Ti_rotated.fasta",format="fasta")
vmatchPattern(consensus_TDNAborder,ti,fixed=FALSE)
