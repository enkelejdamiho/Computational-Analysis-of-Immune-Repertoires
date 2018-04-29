####################################################
## 2018 PEGS Boston
## Immune repertoires with R
## Enkelejda Miho
## ETH Zurich & FHNW Switzerland
## enkelejda.miho@ethz.ch
## enkelejda.miho@fhnw.ch
####################################################

## This file walks you through an introduction in R.
## and analysis of NGS of immune repertoires.
## Thanks to Victor Greiff and Stephen Pettigrew
## and other instructors for letting me pull from 
## their materials.

###################  PART I  #######################
####################################################
##### Quick start in R and immune repertoires ######
####################################################

## Set working directory where you keep the folder
## of the analysis. For example
setwd("~/Desktop/PEGS_2017_EnkelejdaMiho/code")

## Install R packages

source("http://www.bioconductor.org/biocLite.R")
biocLite()

install.packages(c("ggplot2", "knitr"),
                 dependencies = TRUE,
                 repos = "http://cran.us.r-project.org")



## Load R packages

library("ggplot2")
library("knitr")


## ----help_in_R

apropos("mean") # Returns the names of all objects in the search

example("mean") # Examples part of R's online help

help("mean") # Require help regarding the function "mean", equivalent to "?mean". 
# Documentation on a topic with name (typically, an R object or a data set) 
# can be displayed by either help("name") or ?name.

?mean # Access the documentation on a topic with name (e.g. "mean")
?plot # Access the documentation on a topic (e.g. "plot")

??mean # Search the Help System 


## ----load_all_packages

library(pcaMethods)
library(seqinr)
library(RColorBrewer)
library(xtable)
library(plyr)
library(ggplot2)
library(ShortRead)
library(grid)
library(reshape)
library(ape)
library(phylotools)
library(stringr)
library(gridExtra)
library(hexbin)
library(data.table)
library(VennDiagram)
library(scales)
library(fastmatch)
library(HDMD)
library(Biobase)
library(Biostrings)
library(stringdist)
library(ConsensusClusterPlus)
library(Hmisc)
library(gplots)
library(corrgram)
library(igraph)
library(NMF)


## ----data_structures_vector, message=FALSE, warning=FALSE----------------

# Put elements into a numeric vector with the c function, 
# "combine" or "concatenate"
cdr3_lengths <- c(11, 15, 22, 12, 17, 20, 19, 12, 21, 19)

# Store the mean as its own object
mean_cdr3_lengths <- mean(cdr3_lengths) # try sum(), range(), max()

# Print the result
print(mean_cdr3_lengths)

# set.seed makes the sampling reproducible
set.seed(1)

### Character vector of 10 random sequences
cdr3_sequences <- replicate(10,
paste0("CAR", paste(sample(unlist(strsplit('ACDEFGHIKLMNPQRSTVWY', "")),
sample(4:20,1)), collapse = ""), "W"))

# Returns the length of the generated CDR3s
nchar(cdr3_sequences)

# Over-write an object
# cdr3_lengths <- nchar(cdr3_sequences)

cdr3_lengths

nchar(cdr3_sequences)

### Logical vector
cdr3_lengths %in% nchar(cdr3_sequences) 

# NOT: !x
cdr3_lengths[!cdr3_lengths==22] # all values excluding 22

# AND: x&y
# OR: x|y
# XOR, indicates elementwise exclusive OR: xor(x,y)

# Mixed vector and missing data
print(cdr3_sequences[2])

set.seed(1)
mixed_vector <- c(c(rep("A",3), "B"),
                  sample(cdr3_sequences, 4),
                  seq(3, 13, 3))

print(as.numeric(mixed_vector))

# Results in "NA" because not all elements are numeric
mean(as.numeric(mixed_vector)) 

# Remove NA 
mean(as.numeric(mixed_vector), na.rm = T)

# "intersect" performs set intersection
which_lengths_overlap <- intersect(cdr3_lengths, nchar(cdr3_sequences)) 

print(which_lengths_overlap) 


## ----visualization_venn

### Produce figure of overlapping CDR3 lengths in the two repertoires

# install.packages("VennDiagram") # Uncomment to install package "VennDiagram"

library(VennDiagram) # Load package "VennDiagram"

# Output figure 
# pdf("figure/overlap_cdr3_lengths.pdf") # Uncomment this line to produce file
venn_r1_r2 <- venn.diagram(list("Repertoire 1" = cdr3_lengths, 
                                "Repertoire 2" = nchar(cdr3_sequences)), 
                                filename = NULL, print.mode = "raw", 
                                #set print.mode = "perc", for percentage display
                                fill = c("blue",  "orange"), 
                                cex = 1, fontfamily = "sans", 
                                cat.cex = 1, cat.dist = c(0.1,0.1),
                                cat.fontface = "bold", 
                                cat.fontfamily = "sans", margin = 0.1)

grid.draw(venn_r1_r2)
# dev.off() # Uncomment this line to produce the pdf file

# Use ?venn.diagram for a detailed explanation of the arguments used


## ----data_structures_matrix----------------------------------------------

matrix_example <- matrix(mixed_vector, nrow = 4, ncol = 3)
matrix_example
                        
# To order the filling of the matrix by row, set byrow = TRUE
matrix_by_row <- matrix(mixed_vector, nrow = 4, ncol = 3,
                        byrow = T)

matrix_by_row


## ----data_structures_dataframes------------------------------------------

# Transform matrix into dataframe
dataframe_example <- as.data.frame(matrix_example)

# Set row and column names
colnames(dataframe_example) <- c("Sample", "CDR3", "Abundance")
rownames(dataframe_example) <- c(1:nrow(dataframe_example))

dataframe_example

# Add new column
dataframe_example$Vgene <- c("V1-2", "V1-69", "V3-46", "V5-1")


## ----data_structures_arrays----------------------------------------------
array_example <- array(0, dim = c(2, 3, 2))

## ----data_structures_list------------------------------------------------

list_example <- list(cdr3_lengths, cdr3_sequences, dataframe_example,
                matrix_example, array_example)

list_example


## ----accessing_data------------------------------------------------------

cdr3_lengths

### Algebric manipulation
# Algebric operations are vectorized
cdr3_lengths + 1 # Result is each element of the vector +1:

# More operations - ,  /,  %*%

### Extract specific values
cdr3_lengths[3]

cdr3_lengths[2:5] # Sequence of values

matrix_example[1,] # Returns first row

matrix_example[,2] # Returns second column

# Elements in column "CDR3", equivalent to dataframe_example[,"CDR3"]
dataframe_example$CDR3 

colnames(dataframe_example) # Returns column names

names(dataframe_example) # Returns column names

list_example[[2]] # Returns 5th element of list

### Extract unique values
unique(cdr3_lengths)

### Check data: class, length, summary, structure and head
class(cdr3_lengths) # data class (numeric)

length(cdr3_lengths) # vector length

summary(cdr3_lengths) # summary statistics

str(dataframe_example) # Try ?str # Internal structure of an R object

### Change data class
dataframe_example$Abundance <- as.numeric(dataframe_example$Abundance)

dataframe_example$CDR3 <- as.character(dataframe_example$CDR3)

dataframe_example$Vgene <- c("V1-2", "V1-69", "V3-2", "V5-1")

head(dataframe_example) # first elements of an object


## ----visualization_lineplot_density

par(mar=c(4,4,2,.1),cex.lab=.95,cex.axis=.9,mgp=c(2,1,0),tcl=.1)
### Line plot of CDR3 lengths
# Graph CDR3 in function of their length in amino acids
plot(cdr3_lengths, type="o", col="blue", xlab="Antibody clone", 
     ylab="CDR3 length (a.a.)", ylim=c(0,max(nchar(cdr3_sequences))))

# Graph lengths of the CDR3 sequences simulated with red dashed line and square points
lines(nchar(cdr3_sequences), type="o", pch=22, lty=2, col="red")

# Create a title with a red, bold/italic font
title(main="CDR3 length of antibodies in two samples", col.main="black", font.main=1)

### Density plot of CDR3 lengths
plot(density(cdr3_lengths), xlab="CDR3 length (a.a)", main="", col="blue")
lines(density(nchar(cdr3_sequences)), col="red")


## ----visualization_stack_barchart, fig.height=3.4, fig.width=5-----------
### Stacked barchart of a.a. composition
# pdf("figure/stack_barchart_aafreq.pdf") # uncomment this line to save pdf
seqs <- cdr3_sequences
AAs <- alphabetFrequency(AAStringSet(seqs)) ## alphabetFrequency also works on sets
layout(matrix(c(1,2), 1, 2,  byrow = TRUE), width=c(0.7,0.3)) ## allow to place the legend
nAAs = AAs/rowSums(AAs) ## if you want to normalize by total length, use this
barplot(t(nAAs),col=rainbow(ncol(AAs))) ## yields a stacked barplot, one bar per sequence
plot.new()
legend(x="bottom", legend=rev(colnames(AAs)), fill=rev(rainbow(ncol(AAs))),cex=0.7)
# dev.off()

### Piecharts for each of the CDR3 sequences
# for(i in 1:length(cdr3_sequences)){
#   seq = cdr3_sequences[i]
#   AAs = table(strsplit(seq,"", useBytes=TRUE))
#   pie(AAs, col=rainbow(length(AAs)), main="Residue abundance")
# }


## ----visualization_barplot_boxplot

par(mar=c(4,4,.1,.1),cex.lab=.95,cex.axis=.9,mgp=c(2,.7,0),tcl=.1)

barplot(table(cdr3_lengths))

boxplot(cdr3_lengths, nchar(cdr3_sequences))


## ----visualization_vgene

set.seed(7)

vgene <- sample(dataframe_example$Vgene,30, replace = TRUE)

barplot(prop.table(table(vgene))*100, ylab = "Frequency (%)")


## ----visualization_phylogenies

library(ape)

cdr3_cluster <- hclust(stringDist(as.character(cdr3_sequences), method = "levenshtein"))

cdr3_cluster$labels = as.character(cdr3_sequences)

plot(as.phylo(cdr3_cluster))



## ----r_functions

### If statement
# if(logical boolean){function to perform if the boolean is TRUE}
x <- 100000000

if(x > 10000){
  "x is a really big number"
}

if(x < 10000){
  "x isn't that big"
} ## Nothing happens here because x < 10000 evaluates to false

# Further reading: ifelse statement

### For loop
for(i in 1:10){
  print(i)
} # see while loop

# Further reading: while loop

### apply family
apply(matrix_example, 1, nchar) # Count the characters by row (1)

apply(apply(matrix_example, 1, nchar), 2, median) # Median of the previous by column (2)

sapply(1:10, function(x) x+1) # see apply, lapply, tapply


### Writing a function

our_mean_function <- function(x){
  our_sum <- sum(x)
  our_mean <- our_sum / length(x)
  return(our_mean) 
}

our_mean_function(1:10)



###################  PART II #######################
####################################################
### Immune repertoire analysis and visualization ###
####################################################

## ----cdr3_length_distribution

load("clon90")

mean_cdr3_lengths <- tapply(clon90$cdr3_length, clon90$L1, function(x)
  mean(x, na.rm=TRUE))

median_cdr3_lengths <- tapply(clon90$cdr3_length, clon90$L1, function(x)
  median(x, na.rm=TRUE))

pdf("figure/geom_density_cdr3_lengths_clon90.pdf")
ggplot(clon90, aes(x=cdr3_length, color=L1)) + geom_density(size=1.2) +
  scale_color_manual("Donor", values=c("red", "#56B4E9", "darkgreen",
                     "yellow","#999999", "#E69F00"),
                     labels=c("HIV-1 IAVI donor 17 HC", "HIV-1 IAVI donor 17 LC",
                              "Uninfected donor 1 HC", "Uninfected donor 1 LC",
                              "Uninfected donor 2 HC", "Uninfected donor 2 LC")) +
  theme(axis.text.x = element_text(angle=360, size=5)) +
  theme_bw(base_size=12, base_family = "Helvetica") +
  labs (title = "", y="Density", x="CDR3 lengths") +
  theme(strip.background = element_rect(fill="white", linetype="solid", color="white"),
        strip.text=element_text(face="bold", size=rel(1.2), hjust=0.5, vjust=1)) +
  theme(strip.background=element_rect(fill = "white")) +
  theme(strip.background=element_rect(fill = "white")) +
  theme(plot.title=element_text(family="Helvetica", face="bold", size=rel(1.5),
                                hjust=-0.035, vjust=3.5)) +
  theme(axis.text=element_text(family="Helvetica", size=rel(1.2)),
        axis.title=element_text(family="Helvetica", size=rel(1.2)),
        legend.text=element_text(family="Helvetica", size=rel(1.1)),
        legend.title=element_text(family="Helvetica", face="bold", size=rel(1.2))) +
  theme(axis.title.x = element_text(vjust=-0.5, size=rel(1.2)),
        axis.title.y = element_text(vjust=1, size=rel(1.2))) +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line=element_line(size = 0.7, color = "black"),
        text = element_text(size = 14)) +
  theme(plot.margin = unit(c(0, 0.2, 0.5, 0.2),"cm")) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0)) +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line = element_line(size = 0.7, color = "black"),
        text = element_text(size = 12))+ theme(legend.position=c(0.805,0.85)) +
  geom_vline(xintercept=mean_cdr3_lengths, linetype="dashed", size=1,
             color=c("red", "#56B4E9", "darkgreen", "yellow","#999999", "#E69F00"))
dev.off()




##----vgene_frequency_plot

load("working_hiv1_iavi17.RData")

# Germline genes - V gene frequency analysis

vgene_freq <- prop.table(table(working_hiv1_iavi17$vgene_subgroup))

vgene_plot <- ggplot(data.frame(vgene_freq), aes(x=Var1, y=Freq*100)) +
  theme_bw() + geom_bar(stat="identity", width=0.5, fill="dodgerblue4") +
  labs(x="", y="Frequency (%)") +
  theme(axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1)) +
  theme(axis.text.x=element_text(angle=90, size=20, hjust=1, vjust=0.5)) +
  theme(axis.title.y = element_text(colour="black", size=25, vjust=1.5, hjust=0.4),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=20)) +
  theme(plot.margin = unit(c(12, 12, 12, 12),"points"))

pdf("figure/vgene_frequency.pdf", width=12)
vgene_plot
dev.off()




## ----extraction_techniques

# Load data
load("working_hiv1_iavi17.RData")

# As we have ordered our dataframes by CDR3 ranks,
# we can subset the 100 most frequent clones by
top_frequency_cdr3 <- as.character(head(working_hiv1_iavi17[["cdr3"]], 100))

# In order to select all the columns of the top 100 most frequent
# CDR3s in the original dataset:
all_top_cdr3 <- working_hiv1_iavi17[top_frequency_cdr3,]

# Calculate the median mutations in this subset of CDR3s
median(all_top_cdr3$vregion_mutations)
# [1] 21

# FR2 for CDR3 that have V gene IGHV1-2 in the top 100 most frequent CDR3s
data.frame(fr2=all_top_cdr3[all_top_cdr3$vgene_subgroup=="IGHV1-2",]$fr2)

# A faster way is by using the subset function.
# We subset the fr3 in the original dataframe and V genes
# by the clones that have CDR3 longer than 27 a.a.
subset(working_hiv1_iavi17, cdr3_length > 27, select = c(fr3,vgene_subgroup))

# Select for the "CDR3" with all its characteristics (V gene, CDR3 frequency, etc)
working_hiv1_iavi17[working_hiv1_iavi17$cdr3 == "VRDGAYGCSGASCYFGALGNFVYYYYMDV",]

# Find sequence in dataset
which(working_hiv1_iavi17$cdr3 == "VRDGAYGCSGASCYFGALGNFVYYYYMDV")
# [1] 210

# Find similar sequences
cdr3_distance <- stringDist(as.character(working_hiv1_iavi17$cdr3),
                            method = "levenshtein")

# Assign 1 to all similar sequences (max. 1 a.a. different)
# and zero to all sequences that differ by more than 1 a.a.
cdr3_matrix <- ifelse(as.matrix(cdr3_distance) == 1, 1, 0)
colnames(cdr3_matrix) <- as.character(working_hiv1_iavi17$cdr3)
rownames(cdr3_matrix) <- as.character(working_hiv1_iavi17$cdr3)

which(cdr3_matrix[rownames(cdr3_matrix)=="ARRGIAGPDYYSYHGLDV"]==1)
# [1]   52  106  122 1186 1436 1437 1438 1552


## ----save_output

# Example loading presaved data

load("working_hiv1_iavi17.RData")

# Example of saving results as RData

subset_fr3 <- subset(working_hiv1_iavi17, cdr3_length > 27, select = c(fr3,vgene_subgroup))

save(subset_fr3, file="subset_fr3.RData")

# Example of writing results as CSV

top_cdr3 <- as.character(head(working_hiv1_iavi17[["cdr3"]], 100))

write.csv(top_cdr3, "top_cdr3.csv")

# Example of loading and reading-in previous results

load("subset_fr3.RData")

# Print the first row of the data

print(subset_fr3[1,])

nrow(subset_fr3)



## ----visualization

# Load clones based on 90 % CDR3 a.a. identity
load("clon90")

library(e1071)

skewness(clon90[clon90$L1=="HIV-1 IAVI donor 17 HC",]$mean_mutations)
# [1] 1.301486

#placeholder plot - prints nothing at all
empty <- ggplot()+geom_point(aes(1,1), colour="white") +
     theme(
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.border = element_blank(),
       panel.background = element_blank(),
       axis.title.x = element_blank(),
       axis.title.y = element_blank(),
       axis.text.x = element_blank(),
       axis.text.y = element_blank(),
       axis.ticks = element_blank()
     )

#scatterplot of x and y variables
scatter <- ggplot(clon90,aes(cdr3_length, mean_mutations)) +
geom_point(aes(color=L1)) +
scale_color_manual(values = c("red", "#56B4E9", "darkgreen",
                              "yellow","#999999", "#E69F00")) +
theme(legend.position=c(1,1),legend.justification=c(1,1))
  # scatter <- hc_abundance_mutations_plot

#marginal density of x - plot on top
plot_top <- ggplot(clon90,aes(cdr3_length, fill=L1)) +
  geom_density(alpha=.5) +
  scale_fill_manual(values = c("red", "#56B4E9", "darkgreen",
                               "yellow","#999999", "#E69F00")) +
  theme(legend.position = "none")

#marginal density of y - plot on the right
plot_right <- ggplot(clon90,aes(mean_mutations, fill=L1)) +
  geom_density(alpha=.5) +
  coord_flip() +
  scale_fill_manual(values = c("red", "#56B4E9", "darkgreen",
                               "yellow","#999999", "#E69F00")) +
  theme(legend.position = "none")

# Arrange the plots together, with appropriate
# height and width for each row and column
pdf("figure/cdr3_length_mean_mutations.pdf")
grid.arrange(plot_top, empty, scatter, plot_right,
             ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()



## ----networks

load("working_hiv1_iavi17.RData")

library(igraph)

cdr3 <- as.character(working_hiv1_iavi17$cdr3)
cdr3_dist <- stringDist(cdr3, method = "levenshtein")
cdr3_mat <- as.matrix(cdr3_dist)
cdr3_bol <- cdr3_mat
cdr3_bol[cdr3_bol<=1] <- 1
cdr3_bol[cdr3_bol>1] <- 0
colnames(cdr3_bol) <- cdr3
rownames(cdr3_bol) <- cdr3

# CDR3 corresponding to the top aligned VDJ with the bNAb database
as.character(working_hiv1_iavi17[working_hiv1_iavi17$vdj_region ==
"AAMAGTSGPGLVKPSETLSVTCSVSGDSMNNYYWTWIRQSPGKGLEWIGYISDRPSATYNPSLNSR
VVISRDTSKNQLSLKLNSVTPADTAVYYCATARRGQRIYGEVSFGEFFYYYSMDVWGKGTAVTVSS",]$cdr3)

# The entire row
working_hiv1_iavi17[working_hiv1_iavi17$vdj_region==
"AAMAGTSGPGLVKPSETLSVTCSVSGDSMNNYYWTWIRQSPGKGLEWIGYISDRPSATYNPSLNSRVV
ISRDTSKNQLSLKLNSVTPADTAVYYCATARRGQRIYGEVSFGEFFYYYSMDVWGKGTAVTVSS",]

cdr3_graph <- igraph::simplify(graph.adjacency(cdr3_bol, weighted=T, mode = "undirected"))

# Name of the first node of the graph
V(cdr3_graph)$name[1]

nodesize <- ifelse(V(cdr3_graph)$name ==
as.character(working_hiv1_iavi17[working_hiv1_iavi17$vdj_region ==
"AAMAGTSGPGLVKPSETLSVTCSVSGDSMNNYYWTWIRQSPGKGLEWIGYISDRPSATYNPSLNSR
VVISRDTSKNQLSLKLNSVTPADTAVYYCATARRGQRIYGEVSFGEFFYYYSMDVWGKGTAVTVSS",]$cdr3),
3, working_hiv1_iavi17$vregion_mutations/11)

# Mark in red the top hit (best aligned) with bNAb database sequences
nodecolor <- ifelse(V(cdr3_graph)$name ==
as.character(working_hiv1_iavi17[working_hiv1_iavi17$vdj_region ==
"AAMAGTSGPGLVKPSETLSVTCSVSGDSMNNYYWTWIRQSPGKGLEWIGYISDRPSATYNPSLNSR
VVISRDTSKNQLSLKLNSVTPADTAVYYCATARRGQRIYGEVSFGEFFYYYSMDVWGKGTAVTVSS",]$cdr3),
"red", "grey")

pdf("figure/cdr3_network.pdf")
set.seed(11)
plot(cdr3_graph, vertex.frame.color=NA,
     layout=layout_with_fr(cdr3_graph, grid="nogrid", niter=200),
     vertex.label=NA, vertex.color=nodecolor, edge.width=1,
     vertex.size=nodesize, edge.color = "black")
dev.off()

