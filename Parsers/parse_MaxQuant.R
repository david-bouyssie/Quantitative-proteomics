library(stringr)
library(readxl)

# quantif_file = "../Data/Quanti/Raw/MaxQuant/proteinGroups.xlsx"
# output_file = "../Data/Quanti/Parsed_data/Pascaline/proteinGroupshaha.csv"
# source = "MaxQuant"
# intensities_type = "Intensity"
# sheet = 1

# Load data
if(grepl(".xlsx",quantif_file)){
  proteinGroupsInput = read_excel(quantif_file, sheet, col_names = TRUE)
}else{
  proteinGroupsInput <- read.table(quantif_file, header = T, stringsAsFactors = F, sep = "\t")
}

# Extract data
Id = proteinGroupsInput[,grepl("Protein.IDs",names(proteinGroupsInput))]
colnames(Id) = "Id"
Accession = Id
colnames(Accession) = "Accession"
Gene_name = proteinGroupsInput[,grepl("Gene.names",names(proteinGroupsInput))]
if(length(Gene_name)>0){
  colnames(Gene_name) = "Gene_name"
}

intensities = proteinGroupsInput[,grepl(paste0(intensities_type,".+"),names(proteinGroupsInput))]
if(intensities_type=="LFQ"){
  colnames(intensities) = paste0("Intensity_",sub(paste0(".*",intensities_type,".intensity."), "", colnames(intensities)))
}else{
  colnames(intensities) = paste0("Intensity_",sub(paste0(".*",intensities_type,"."), "", colnames(intensities)))
}
intensities[intensities == 0] <- NA

identification_types = proteinGroupsInput[,grepl("Identification.type.",names(proteinGroupsInput))]
identification_types[identification_types == ""] <- NA
colnames(identification_types) = paste0("Identification_type_",sub("Identification.type.","",colnames(identification_types)))

# Specific peptides
specific_peptides = proteinGroupsInput[grepl("^Unique.peptides$",names(proteinGroupsInput))]
colnames(specific_peptides)="Specific_peptides"

# Build dataframe
proteinGroupsOutput = cbind(Id,Accession,Gene_name,identification_types,specific_peptides,intensities)
write.table(proteinGroupsOutput,output_file,row.names=FALSE,sep="\t")
