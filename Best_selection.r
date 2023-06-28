###article title: What’s left in the tank? Identification of non-ascribed aquarium’s coral collections with DNA barcodes as part of an integrated diagnostic approach
###journal name: Conservation Genetics Resources
###author names: Luigi Colin; Daniel Abed-Navandi; Dalia A. Conde; Jamie Craggs; Ana Rita da Silva; Max Janse; Björn Källström; Alexander Pearce-Kelly; Chris Yesson
###Corresponding author: luigi.colin@ioz.ac.uk
###
###Confidence Values and species two region preference selection
###Based on quantile distribution of references sequence used for Barcoding gap analysis 
###mtCR: <98, "Diff sp"; 98 < X < 99.5, "<50%"; 99.5 < X < 100, "50-99%"; X = 100, ">99%"
###PaxC: <90, "Diff sp"; 90 < X < 98.5, "<50%"; 98.5 < X < 99.998, "50-99%"; X > 99.998, ">99%"

library("xlsx") # Only if clean file is xlsx

mtCR <- read.xlsx2("~/Analysis/ZSL/Barcoding/Blast/mtCR.xlsx", sheetIndex = 1) # Read id data from previous "blast match" code
PaxC <- read.xlsx2("~/Analysis/ZSL/Barcoding/Blast/PaxC.xlsx", sheetIndex = 1) # Read id data from previous "blast match" code

mtCR$Confidence <- 0
PaxC$Confidence <- 0

mtCR <- mtCR[!grepl(pattern = "NA", mtCR$Perc.Ident),]

for (i in 1:length(mtCR$Perc.Ident)) {
  x <- as.numeric(as.character(mtCR$Perc.Ident[i]))
  ### mtCR: <98, "Diff sp"; 98 < X < 99.5, "<50%"; 99.5 < X < 100, "50-99%"; X = 100, ">99%"
  if (x < 98) {
    mtCR$Confidence[i] <- "Diff sp."
  }
  if ((x >= 98) && (x < 99.5)) {
    mtCR$Confidence[i] <- "<50%"
  }
  if (x >= 99.5) {
    mtCR$Confidence[i] <- "50-99%"
  }
  if (x == 100) {
    mtCR$Confidence[i] <- ">99%"
  }
}

# Remove "NA"
PaxC <- PaxC[!grepl(pattern = "NA", PaxC$Perc.Ident),]

for (i in 1:length(PaxC$Perc.Ident)) {
  x <- as.numeric(as.character(PaxC$Perc.Ident[i]))
  ### PaxC: <90, "Diff sp"; 90 < X < 98.5, "<50%"; 98.5 < X < 99.998, "50-99%"; X > 99.998, ">99%"
  if (x < 90) {
    PaxC$Confidence[i] <- "Diff sp."
  }
  if ((x >= 90) && (x < 98.5)) {
    PaxC$Confidence[i] <- "<50%"
  }
  if ((x >= 98.5) && (x < 99.998)) {
    PaxC$Confidence[i] <- "50-99%"
  }
  if (x >= 99.998) {
    PaxC$Confidence[i] <- ">99%"
  }
}

id <- merge(mtCR, PaxC, by = "QueryID", suffixes = c("-mtCR", "-PaxC"), all.x = TRUE, all.y = TRUE, sort = TRUE)
id2 <- merge(mtCR, PaxC, by = "QueryID", suffixes = c("-mtCR", "-PaxC"), sort = TRUE)

id.sim <- id[, c("QueryID", "Organism-mtCR", "Organism-PaxC", "GR-mtCR", "GR-PaxC", "Confidence-mtCR", "Confidence-PaxC", "Perc.Ident-mtCR", "Perc.Ident-PaxC", "Bits-mtCR", "Bits-PaxC")]
id.full <- id[, c("QueryID", "Organism-mtCR", "Organism-PaxC", "GR-mtCR", "GR-PaxC", "Confidence-mtCR", "Confidence-PaxC", "SubjectID-mtCR", "SubjectID-PaxC", "Perc.Ident-mtCR", "Perc.Ident-PaxC", "Bits-mtCR", "Bits-PaxC", "Alignment.Length-mtCR", "Alignment.Length-PaxC", "Mismatches-mtCR", "Mismatches-PaxC")]

# Write
write.table(id.sim, "~/Analysis/ZSL/Barcoding/Blast/Temp/ID-simple.csv", sep = ",", row.names = FALSE)
write.table(id.full, "~/Analysis/ZSL/Barcoding/Blast/Temp/ID-full.csv", sep = ",", row.names = FALSE)
write.table(id2, "~/Analysis/ZSL/Barcoding/Blast/Temp/Overlap-full.csv", sep = ",", row.names = FALSE)

id2 <- id.full

id2$ID.test <- "not done"
id2$notes.test <- ""

# Conflict solving
### Choose highest confidence where they don't match at all
## Remove NA values
id2[is.na(id2)] <- 0

# Define confidence classes as numerical
class <- data.frame(Class = unique(id2$`Confidence-mtCR`), Value = c(4, 3, 2, 1, 0))
id2$CCV.mtCR <- 0
id2$CCV.PaxC <- 0

# Assign class values to both genes
for (i in 1:length(id2$QueryID)) {
  id2$CCV.mtCR[i] <- class$Value[class$Class == id2$`Confidence-mtCR`[i]]
  id2$CCV.PaxC[i] <- class$Value[class$Class == id2$`Confidence-PaxC`[i]]
}

# Choose
for (i in 1:length(id2$QueryID)) {
  if (id2$`Confidence-mtCR`[i] == id2$`Confidence-PaxC`[i]) {
    id2$ID.test[i] <- "no choice"
    id2$notes.test[i] <- "Same confidence, conflict not resolved"
  }
  if (id2$CCV.mtCR[i] > id2$CCV.PaxC[i]) {
    id2$ID.test[i] <- as.vector(id2$`Organism-mtCR`[i])
    id2$notes.test[i] <- "mtCR higher confidence"
  }
  if (id2$CCV.mtCR[i] < id2$CCV.PaxC[i]) {
    id2$ID.test[i] <- as.vector(id2$`Organism-PaxC`[i])
    id2$notes.test[i] <- "PaxC higher confidence"
  }
}
#----confidence ----

# If mtCR confidence is "Diff sp.", choose PaxC ID
for (i in 1:length(id2$QueryID)) {
  if (id2$`Confidence-mtCR`[i] == "Diff sp.") {
    id2$ID.test[i] <- as.vector(id2$`Organism-PaxC`[i])
    id2$notes.test[i] <- "No confidence in mtCR data, best ID PaxC"
  }
}

# If PaxC confidence is "Diff sp.", choose mtCR ID
for (i in 1:length(id2$QueryID)) {
  if (id2$`Confidence-PaxC`[i] == "Diff sp.") {
    id2$ID.test[i] <- as.vector(id2$`Organism-mtCR`[i])
    id2$notes.test[i] <- "No confidence in PaxC data, best ID mtCR"
  }
}

#----match gr ----

# Confirm one PaxC ID with mtCR data at the group species level
for (i in 1:length(id2$QueryID)) {
  if (is.na(id$`GR-mtCR`[i])) {
    next()
  }
  if (is.na(id$`GR-PaxC`[i])) {
    next()
  }
  if (grepl(pattern = "not determined", id2$`GR-mtCR`[i])) {
    next()
  }
  
  patt <- as.vector(id2$`GR-PaxC`[i])
  if (grepl(pattern = patt, id2$`GR-mtCR`[i])) {
    id2$ID.test[i] <- as.vector(patt)
    id2$notes.test[i] <- "Confirmation of one PaxC ID with mtCR data at group sp. level"
  }
}

# Confirm one mtCR ID with PaxC data at the group species level
for (i in 1:length(id2$QueryID)) {
  if (is.na(id$`GR-mtCR`[i])) {
    next()
  }
  if (is.na(id$`GR-PaxC`[i])) {
    next()
  }
  patt <- as.vector(id2$`GR-mtCR`[i])
  if (grepl(pattern = "not determined", id2$`GR-PaxC`[i])) {
    next()
  }
  if (grepl(pattern = patt, id2$`GR-PaxC`[i])) {
    id2$ID.test[i] <- as.vector(patt)
    id2$notes.test[i] <- "Confirmation of one mtCR ID with PaxC data at group sp. level"
  }
}

##----ID matching sp gr----

# Match at the group species level
for (i in 1:length(id2$QueryID)) {
  if (is.na(id$`GR-mtCR`[i])) {
    next()
  }
  if (is.na(id$`GR-PaxC`[i])) {
    next()
  }
  
  if (as.vector(id2$`GR-mtCR`[i]) == as.vector(id2$`GR-PaxC`[i])) {
    id2$ID.test[i] <- as.vector(id2$`GR-mtCR`[i])
    id2$notes.test[i] <- "Match at group sp. level"
  }
}

# Confirm one mtCR ID with PaxC data
for (i in 1:length(id2$QueryID)) {
  if (is.na(id$`Organism-mtCR`[i])) {
    next()
  }
  if (is.na(id$`Organism-PaxC`[i])) {
    next()
  }
  patt <- as.vector(id2$`Organism-PaxC`[i])
  if (grepl(pattern = patt, id2$`Organism-mtCR`[i])) {
    id2$ID.test[i] <- as.vector(patt)
    id2$notes.test[i] <- "Confirmation of one mtCR ID with PaxC data"
    id2$Nloop[i] <- id2$Nloop[i] + 1
  }
}

# Confirm one PaxC ID with mtCR data
for (i in 1:length(id2$QueryID)) {
  if (is.na(id$`Organism-mtCR`[i])) {
    next()
  }
  if (is.na(id$`Organism-PaxC`[i])) {
    next()
  }
  patt <- as.vector(id2$`Organism-mtCR`[i])
  if (grepl(pattern = patt, id2$`Organism-PaxC`[i])) {
    id2$ID.test[i] <- as.vector(patt)
    id2$notes.test[i] <- "Confirmation of one PaxC ID with mtCR data"
    id2$Nloop[i] <- id2$Nloop[i] + 1
  }
}

##-----ID matching sp-----

# Match of species ID
for (i in 1:length(id2$QueryID)) {
  if (is.na(id$`Organism-mtCR`[i])) {
    next()
  }
  if (is.na(id$`Organism-PaxC`[i])) {
    next()
  }
  if (as.vector(id2$`Organism-mtCR`[i]) == as.vector(id2$`Organism-PaxC`[i])) {
    id2$ID.test[i] <- as.vector(id2$`Organism-mtCR`[i])
    id2$notes.test[i] <- "Match of species ID"
  }
}

# No mtCR data, choose PaxC ID
for (i in 1:length(id2$QueryID)) {
  if (id2$`Confidence-mtCR`[i] == 0) {
    id2$ID.test[i] <- as.vector(id2$`Organism-PaxC`[i])
    id2$notes.test[i] <- "no mtCR data, best ID PaxC"
  }
  if ((id2$`Confidence-mtCR`[i] == 0) && (id2$CCV.PaxC[i] == 1)) {
    id2$ID.test[i] <- as.vector(id2$`Organism-PaxC`[i])
    id2$notes.test[i] <- "no mtCR data, no confidence in PaxC ID"
  }
  if (id2$`Confidence-PaxC`[i] == 0) {
    id2$ID.test[i] <- as.vector(id2$`Organism-mtCR`[i])
    id2$notes.test[i] <- "no PaxC data, best ID mtCR"
  }
  if ((id2$`Confidence-PaxC`[i] == 0) && (id2$CCV.mtCR[i] == 1)) {
    id2$ID.test[i] <- as.vector(id2$`Organism-mtCR`[i])
    id2$notes.test[i] <- "no PaxC data, no confidence in mtCR ID"
  }
}

###-----write file----

# Write the updated table to a CSV file
# write.table(id2, "~/Analysis/ZSL/Barcoding/Blast/Temp/ID-table.csv", sep = ",", row.names = F)

# Write the updated table to an Excel file
write.xlsx2(id2, "~/Analysis/ZSL/Barcoding/Blast/Temp/ID-table.xlsx", row.names = F)