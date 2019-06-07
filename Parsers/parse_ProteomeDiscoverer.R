library(stringr)
library(readxl)

# Load data
if(grepl(".xlsx",quantif_file)){
  proteinGroupsInput = read_excel(quantif_file, 1, col_names = TRUE)
}else{
  proteinGroupsInput <- read.table(quantif_file, header = T, stringsAsFactors = F, sep = "\t")
}

# Extract data
Id = proteinGroupsInput[,grepl("Protein.Group.IDs",names(proteinGroupsInput))]
Accession = proteinGroupsInput[,grepl("^Accession$",names(proteinGroupsInput))]

intensities = proteinGroupsInput[,grepl(paste0(intensities_type," .+Sample"),names(proteinGroupsInput))]
colnames(intensities) = paste0("Intensity_",sub(".*Sample ", "", colnames(intensities)))

identification_types = proteinGroupsInput[,grepl("Found in Sample in",names(proteinGroupsInput))]
identification_types[identification_types=="Not Found"] <- NA
identification_types[identification_types=="High"] <- "By matching"
identification_types[identification_types=="Peak Found"] <- "By MS/MS"
colnames(identification_types) = paste0("Identification_type_",sub("Found in Sample in .* Sample ","",colnames(identification_types)))

# Specific peptides
specific_peptides = proteinGroupsInput[grepl("Number.of.Unique.Peptides",names(proteinGroupsInput))]
colnames(specific_peptides)="Specific_peptides"

# Build dataframe
proteinGroupsOutput = cbind(Id,Accession,identification_types,specific_peptides,intensities)
write.csv(proteinGroupsOutput,output_file,row.names=FALSE)
