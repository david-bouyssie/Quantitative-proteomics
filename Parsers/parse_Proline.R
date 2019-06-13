library(stringr)
library(readxl)

# Parameters set in launcher script (see "Lauchers" section)
# quantif_file = quantification file to parse
# output_file = output table
# intensities_type = abundance/raw_abundance
# sheet = sheet number if quantif_file is a xlsx file

# Load data
if(grepl(".xlsx",quantif_file)){
  proteinGroupsInput = read_excel(quantif_file, sheet, col_names = TRUE)
}else{
  proteinGroupsInput <- read.table(quantif_file, header = T, stringsAsFactors = F, sep = "\t")
}


# Extract data
Id = proteinGroupsInput[grepl("id$",names(proteinGroupsInput))][1]
colnames(Id)="Id"


Accession = proteinGroupsInput[grepl("accession$",names(proteinGroupsInput))][1]
colnames(Accession)="Accession"

Gene_name = proteinGroupsInput[grepl("gene_name",names(proteinGroupsInput))]
if(length(Gene_name)>0){
  colnames(Gene_name)="Gene_name"
}

intensities = proteinGroupsInput[grepl(paste0("^",intensities_type,".+"),names(proteinGroupsInput))]
colnames(intensities) = paste0("Intensity",sub(intensities_type, "", colnames(intensities)))

identification_types = proteinGroupsInput[grepl("psm_count",names(proteinGroupsInput))]
samples = sub("psm_count_", "", colnames(identification_types))
for(sample in samples){
  identification_type = identification_types[grepl(paste0("_",sample),names(identification_types))]
  colname = colnames(identification_type)[1]
  identification_type[[colname]][identification_type[[colname]]>0 & !(is.na(identification_type[[colname]]))] <- "By MS/MS"
  identification_type[[colname]][identification_type[[colname]]==0 & !(is.na(identification_type[[colname]]))] <- "By matching"
  identification_type[[colname]][is.na(identification_type[[colname]])] <- NA
  identification_types[grepl(paste0("_",sample),names(identification_types))] = identification_type[[colname]]
}
colnames(identification_types) = paste0("Identification_type_",samples)

specific_peptides = proteinGroupsInput[grepl("specific_peptide_matches",names(proteinGroupsInput))]
if(length(specific_peptides)>0){
  colnames(specific_peptides)="Specific_peptides"
}

# Build dataframe
proteinGroupsOutput = cbind(Id,Accession,Gene_name,identification_types,specific_peptides,intensities)
write.table(proteinGroupsOutput,output_file,row.names=FALSE,sep="\t")
