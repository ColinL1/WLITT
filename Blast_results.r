### article title: What's left in the tank? Identification of non-ascribed aquarium's coral collections with DNA barcodes as part of an integrated diagnostic approach
### journal name:Conservation Genetics Resources
### author names: Luigi Colin; Daniel Abed-Navandi; Dalia A. Conde; Jamie Craggs; Ana Rita da Silva; Max Janse; Bjorn Kallstrom; Alexander Pearce-Kelly; Chris yesson
### Corresponding author: luigi.colin@ioz.ac.uk

### Blast match and "most probable" ID selection

library("rBLAST")
require("rentrez") # sp name
library("xlsx") # to read sp.group file

setwd("~/Analysis/") # set working directory

# ---- Set up ----
## define region name to be used in file names
gene <- "PaxC" # "rand", "PaxC", "mtCR"
## suffix for csv output file
suff <- "" # "-gr"; "-test"
## define folder path / file names
## parameter filepath
f.path <- c("ZSL/Barcoding/Blast/", "ZSL/Barcoding/Blast/Multi-ID/", "ZSL/Barcoding/Blast/Temp/") # first where to save main result, second where to save multi ID data, third where to save temp file (to debug/check final output)
## parameter namestem
f.result <- paste(gene, suff, ".csv", sep = "") # change name of output file; don't change .csv extension without changing the rest of the code.
f.mm <- c(paste(gene, "-bit", suff, ".csv", sep = ""), paste(gene, "-PID", suff, ".csv", sep = ""), paste(gene, "-misc", suff, ".csv", sep = "")) # temp file names Bit == highest Bit score match, PID == Highest % identity match, MM == Lowest number of mismatches; misc == general information to debug (rarely used)

## Load gr.sp file; only necessary if multi ID are present (very likely); if not, change code to exclude group ID
gr.sp <- read.xlsx("~/Analysis/ZSL/Barcoding/Taxonomy/Group-sp_WQM.xlsx", sheetIndex = 1) # data from "Revision and catalogue of worldwide staghorn corals Acropora and Isopora (Scleractinia: Acroporidae) in the Museum of Tropical Queensland"
gr.sp <- gr.sp[, c(3, 1)] # reorder/clean dataframe

## create blast database if not already made; Fasta file downloaded (Nov 15, 2019) from https://www.ncbi.nlm.nih.gov/; Text query: '"Acroporidae"[Organism]'.
# makeblastdb("~/Analysis/ZSL/All-acroporidea.fasta", dbtype = "nucl", args = "max_file_sz 4000000000")

# Load db (no extension)
bl <- blast(db = "~/Analysis/ZSL/Barcoding/Adb/ADB")
bl

## load sequences to ID (non-aligned fasta file)
seq <- readDNAStringSet("~/Analysis/ZSL/Barcoding/Acropora/Sequences/PaxC-Plate1.fasta", format = "fasta")

## set proportion of minimum length of alignment match (i.e.: 2/3 of seq query length for mtCR)
mlm <- 2/3

# ---- Loop Blast match and choices (for single gene) ----
# enterez query limit to 3 per second. If name sp from local record, then remove Sys.sleep(0.4)
for (j in 1:length(seq)) {
  print(paste("number", j, "out of", length(seq), "started", sep = " "))
  ## Blast match of seq J against references 
  cl <- predict(bl, seq[j,])
  cl <- cl[which(cl$Alignment.Length >= (width(seq[j, ]) * mlm)), ] # remove matches with length less than mlm (set above) of the seq j
  ### skip if object cl empty (no match above minimum set above) + write on csv "skipped j iteration"
  if (dim(cl)[1] == 0) {
    nt.mtc <- data.frame()
    nt.mtc <- data.frame(QueryID = names(seq[j]), SubjectID = "skipped for too short length match", Organism = "NA", Max.Bits = "NA", max.Perc.Ident = "NA", low.Mismatches = "NA", Alignment.Length = "NA") # nice to add length of query here.
    write.table(nt.mtc, file = paste(f.path[1], f.result, sep = ""), append = TRUE, row.names = FALSE, sep = ",")
    print(paste("number", j, "out of", length(seq), "skipped", sep = " "))
    next()
  }
  
  cl2 <- cl[order(cl$Bits, cl$Perc.Ident, decreasing = TRUE), ] # order results by highest Bit  
  # Divide result into highest Bit score and Highest %identity 
  a <- cl2[which(cl2$Perc.Ident == max(cl2$Perc.Iden)), ] # %identity
  b <- cl2[which(cl2$Bits == max(cl2$Bits)), ] # Bit score

  # check if a, b have differences ==> if not, do only one (less server requests ==> less time, fewer failures)
  if (!identical(a, b)) {
    print("a != b")
    # add sp names field (to be filled)
    a$Organism <- as.character(0)
    b$Organism <- as.character(0)

    ## entrez loops add organism name to accession number
    ## highest %identity ## fill in a
    print(paste("adding species names to a (highest Perc.Ident) seq n=", j, sep = "")) # for reference 
    for (i in 1:length(a$SubjectID)) {
      print(as.vector(a$SubjectID[i]))
      # construct a string of your query using standard search terms
      myquery <- paste(a$SubjectID[i])
      # do a search on entrez/genbank nucleotide database
      n <- entrez_search(db = "nucleotide", term = myquery, retmax = 100)
      # get species name
      t <- entrez_summary(db = "nucleotide", id = n$ids, retmax = 100)

      # fill column with number of hits
      a$Organism[i] <- t$organism
      Sys.sleep(0.4)
    }
    ## Highest bitscore ## fill in b
    print(paste("adding species names to b (Highest Bit score) seq n=", j, sep = "")) # for reference 
    for (i in 1:length(b$SubjectID)) {
      print(as.vector(b$SubjectID[i]))
      # construct a string of your query using standard search terms
      myquery <- paste('"', b$SubjectID[i])
      # do a search on entrez/genbank nucleotide database
      n <- entrez_search(db = "nucleotide", term = myquery, retmax = 100)
      # get species name
      t <- entrez_summary(db = "nucleotide", id = n$ids, retmax = 100)
      # fill column with number of hits
      b$Organism[i] <- t$organism
      Sys.sleep(0.4)
    }
  } else {
    if (identical(a, b)) {
      # add sp names field (for references)
      print("a = b")
      a$Organism <- as.character(0)
      ## highest %identity ## fill in a
      print(paste("adding species names to a (highest Perc.Ident) seq n=", j, sep = "")) # for reference 
      for (i in 1:length(a$SubjectID)) {
        print(as.vector(a$SubjectID[i]))
        # construct a string of your query using standard search terms
        myquery <- paste(a$SubjectID[i])
        # do a search on entrez/genbank nucleotide database
        n <- entrez_search(db = "nucleotide", term = myquery, retmax = 100)
        # get species name
        t <- entrez_summary(db = "nucleotide", id = n$ids, retmax = 100)

        # fill column with number of hits
        a$Organism[i] <- t$organism
        Sys.sleep(0.4)
      }
      b <- a # make b from a (ok as they were identical at the beginning)
    }
  }

 ## Order columns a and b
p.id <- a[, c("QueryID", "SubjectID", "Organism", "Bits", "Perc.Ident", "Mismatches", "Alignment.Length")]
bits <- b[, c("QueryID", "SubjectID", "Organism", "Bits", "Perc.Ident", "Mismatches", "Alignment.Length")]

# Cleaning step ==> remove misc name sp. results
p.id <- p.id[!grepl("sp.", p.id$Organism), ] # remove misc name sp. from Highest %identity
bits <- bits[!grepl("sp.", bits$Organism), ] # remove misc name sp. from Highest Bit score

## Check if empty
if ((dim(p.id)[1] == 0) && (dim(bits)[1] != 0)) {
  p.id <- bits # if %identity empty, then copy from Bit score to avoid code errors.
}
if ((dim(bits)[1] == 0) && (dim(p.id)[1] != 0)) {
  bits <- p.id # if Bit score empty, then copy from %identity to avoid code errors.
}
if ((dim(p.id)[1] == 0) && (dim(bits)[1] == 0)) {
  ## If both empty, then save as no-results/match to misc sequence.
  nt.mtc <- data.frame(QueryID = names(seq[j]), SubjectID = "match only to no name sequence", Organism = "NA", Max.Bits = "NA", max.Perc.Ident = "NA", low.Mismatches = "NA", Alignment.Length = "NA") # nice to add length of query here.
  write.table(nt.mtc, file = paste(f.path[1], f.result, sep = ""), append = TRUE, row.names = FALSE, sep = ",")
  print(paste("number", j, "out of", length(seq), "skipped - no name match", sep = " "))
  next()
}

## Add sp-groups to preliminary match
## Integration sp.group
bits$GR <- as.character(0) # add sp.group field
## Fill group species info from loaded file gr.sp
for (i in 1:length(bits$QueryID)) {
  g.s <- as.vector(gr.sp$Species.group[(gr.sp$Species %in% bits[i, ]$Organism)]) # find group matching species
  if (length(g.s) == 1) {
    bits$GR[i] <- g.s # if match, add GR to bits object
  }
  if (length(g.s) == 0) {
    bits$GR[i] <- as.character(paste("not determined", sep = "")) # if no match, add "not determined" to bits object
  }
}

p.id$GR <- as.character(0) # add sp.group field, same as above but for %identity
## Fill group species info from loaded file gr.sp
for (i in 1:length(p.id$QueryID)) {
  g.s <- as.vector(gr.sp$Species.group[(gr.sp$Species %in% p.id[i, ]$Organism)]) # find group matching species
  if (length(g.s) == 1) {
    p.id$GR[i] <- g.s # if match, add GR to p.id object
  }
  if (length(g.s) == 0) {
    p.id$GR[i] <- as.character(paste("not determined", sep = "")) # if no match, add "not determined" to p.id object
  }
}

## Save in /temp for debug/deeper analysis
write.table(p.id, file = paste(f.path[3], f.mm[2], sep = ""), sep = ",", row.names = FALSE, append = TRUE)
## Save in /temp for debug/deeper analysis
write.table(bits, file = paste(f.path[3], f.mm[1], sep = ""), sep = ",", row.names = FALSE, append = TRUE)

##----choose best match ----
##fisrt save result if there is a consensus
if ((length(unique(p.id$Organism)) == 1) && (length(unique(bits$Organism)) == 1) && (p.id[1,]$Organism == bits[1,]$Organism)) {
  write.table(p.id[1,], file = paste(f.path[1], f.result, sep = ""), append = TRUE, row.names = FALSE, sep = ",")
} else { 
  ##save result if 100% match to one species ##uses which.max to select only the 100% match with the highest bit score.
  if (((p.id[1,]$Perc.Ident) == 100) && (length(unique((p.id[grepl((p.id[which.max(p.id$Bits),]$Bits), p.id$Bits),]$Organism))) == 1)) {
    write.table(p.id[1,], file = paste(f.path[1], f.result, sep = ""), append = TRUE, row.names = FALSE, sep = ",")
  } else {
    ##save result if 100% match to unique species group rather than sp. ##prioritize perfect match (100% to group rather than lower % to single species.)
    if ((length(unique((p.id[grepl((p.id[which.max(p.id$Bits),]$Bits), p.id$Bits),]$GR))) == 1) && ((p.id[1,]$Perc.Ident) == 100)) {
      g.id <- data.frame(QueryID = p.id[1,]$QueryID, SubjectID = p.id[1,]$SubjectID, Organism = paste(length(unique(p.id$Organism)), " species (", toString(unique(p.id$Organism)), ")", sep = ""), Bits = p.id[1,]$Bits, Perc.Ident = p.id[1,]$Perc.Ident, Mismatches = p.id[1,]$Mismatches, Alignment.Length = p.id[1,]$Alignment.Length, GR = paste("gr. ", p.id[1,]$GR, sep = ""))
      ##make dataframe with species group "organism" and 100% match info   
      write.table(g.id, file = paste(f.path[1], f.result, sep = ""), append = TRUE, row.names = FALSE, sep = ",")
    } else {
      ##if no 100% match (as above) save result if highest Bit score is a match to one species ##uses which.max to select only the highest %identity.
      if (length(unique(bits[grepl((bits[which.max(bits$Perc.Ident),]$Bits), bits$Bits),]$Organism)) == 1) {
        write.table(bits[1,], file = paste(f.path[1], f.result, sep = ""), append = TRUE, row.names = FALSE, sep = ",")
      } else {
        ## same as above applied to highest Bit score
        if (length(unique(bits[grepl((bits[which.max(bits$Perc.Ident),]$Bits), bits$Bits),]$GR)) == 1) {
          m.id <- data.frame()
          m.id <- data.frame(QueryID = bits[1,]$QueryID, SubjectID = bits[1,]$SubjectID, Organism = paste(length(unique(bits$Organism)), " species (", toString(unique(bits$Organism)), ")", sep = ""), Bits = bits[1,]$Bits, Perc.Ident = bits[1,]$Perc.Ident, Mismatches = bits[1,]$Mismatches, Alignment.Length = bits[1,]$Alignment.Length, GR = paste("gr. ", bits[1,]$GR, sep = ""))
          ##fill dataframe with species group "organism" and best Bit score match info   
          write.table(m.id, file = paste(f.path[1], f.result, sep = ""), append = TRUE, row.names = FALSE, sep = ",")
        } else {
          ##if still no unique species name match ==> save number of group
          if (((p.id[1,]$Perc.Ident) == 100) && (length(unique((p.id[grepl((p.id[which.max(p.id$Bits),]$Bits), p.id$Bits),]$GR))) != 1)) {
            m.id <- data.frame()
            m.id <- data.frame(QueryID = p.id[1,]$QueryID, SubjectID = p.id[1,]$SubjectID, Organism = paste(length(unique(p.id$Organism)), " species (", toString(unique(p.id$Organism)), ")", sep = ""), Bits = p.id[1,]$Bits, Perc.Ident = p.id[1,]$Perc.Ident, Mismatches = p.id[1,]$Mismatches, Alignment.Length = p.id[1,]$Alignment.Length, GR = paste(length(unique(p.id$GR)), " species groups (", toString(unique(p.id$GR)), ")", sep = ""))
            ##fill dataframe with species group number and Highest %identity match info
            write.table(m.id, file = paste(f.path[1], f.result, sep = ""), append = TRUE, row.names = FALSE, sep = ",")
          } else {
            if (((p.id[1,]$Perc.Ident) != 100) && (length(unique(bits[grepl((bits[which.max(bits$Perc.Ident),]$Bits), bits$Bits),]$GR)) != 1)) {
              m.id <- data.frame()
              m.id <- data.frame(QueryID = bits[1,]$QueryID, SubjectID = bits[1,]$SubjectID, Organism = paste(length(unique(bits$Organism)), " species (", toString(unique(bits$Organism)), ")", sep = ""), Bits = bits[1,]$Bits, Perc.Ident = bits[1,]$Perc.Ident, Mismatches = bits[1,]$Mismatches, Alignment.Length = bits[1,]$Alignment.Length, GR = paste(length(unique(bits$GR)), " species groups (", toString(unique(bits$GR)), ")", sep = ""))
              ##fill dataframe with species group number and best Bit score match info
              write.table(m.id, file = paste(f.path[1], f.result, sep = ""), append = TRUE, row.names = FALSE, sep = ",")
            }
          }
          ##save temp file with multiID info
          if (identical(p.id, bits)) {
            write.table(bits, file = paste(f.path[2], f.mm[1], sep = ""), append = TRUE, row.names = FALSE, sep = ",")
          } ##if identical write only Bit score file
          if (!identical(p.id, bits)) {
            ##if not identical write both Bits and p.id file
            write.table(bits, file = paste(f.path[2], f.mm[1], sep = ""), append = TRUE, row.names = FALSE, sep = ",")
            write.table(p.id, file = paste(f.path[2], f.mm[2], sep = ""), append = TRUE, row.names = FALSE, sep = ",")
          }
        }
      }
    }
  }
}

print(paste("number", j, "out of", length(seq), "done!", sep = " "))
}

#---- additional output cleaning for easier reading CSV---- 
##possible to integrate in above script ##easily done out of R.
blast.r <- read.csv(file = paste(f.path[1], f.result, sep = ""))  ##read result file
blast.r <- blast.r[!grepl("QueryID", blast.r$QueryID), ]  #remove extra header
#blast.r <- blast.r[!grepl("skipped", blast.r$SubjectID),] #remove skipped
write.table(blast.r, file = paste(f.path[1], "clean-", f.result, sep = ""), append = FALSE, row.names = FALSE, sep = ",")  ##write result file cleaned.