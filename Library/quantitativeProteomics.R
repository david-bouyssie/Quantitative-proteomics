
#################
# Normalization #
#################

# This function normalize intensities using the ???probabilistic quotient normalization??? method (Dieterle et al, 2006).
normalize <- function(df, reference_sample, method = "pqn"){
  
  # Inputs :
  # df : quantification matrix whom 1rst column contains proteins id
  # reference_sample : name of the column corresponding to the reference sample
  # Output :
  # dataframe containing the ratio values for each sample, and the normalized quantification matrix
  
  # Initialize results dataframe
  ratio_names = paste0(colnames(df)[!colnames(df)%in%colnames(df)[1] & !colnames(df)%in%reference_sample],"_ratio_",reference_sample)
  norm_names = paste0(colnames(df)[-1],"_norm")
  normalized_data = data.frame(matrix(vector(),
                                      nrow(df),
                                      (length(ratio_names)+length(norm_names)+1),
                                      dimnames = list(c(), c("Id",ratio_names,norm_names) )),
                               check.names=FALSE)
  
  # Compute ratios
  normalized_data$Id <- df[[1]]
  reference_intensities = as.numeric(as.vector(unlist(df[as.character(reference_sample)])))
  for (sample in colnames(df)[-1]){
    if(sample!=reference_sample){
      normalized_data[paste0(sample,"_ratio_",reference_sample)]<-as.numeric(as.vector(unlist(df[as.character(sample)])))/reference_intensities
    }
  }
  
  if (method == "pqn"){
    # Compute median of ratios
    median_ratios=c()
    ratios = normalized_data[grepl("_ratio_",colnames(normalized_data))]
    for (ratio in ratios){
      median = median(ratio,na.rm=T)
      median_ratios = c(median_ratios,median)
    }
    
    # Median removal
    i=1
    for (sample in colnames(df)[-1]){
      if(sample==reference_sample){
        normalized_data[paste0(reference_sample,"_norm")]<-reference_intensities
      }else{
        normalized_data[paste0(sample,"_norm")]<-as.numeric(as.vector(unlist(df[as.character(sample)])))/median_ratios[i]
        i=i+1
      }
    }
  }else{
    # Compute median of ratios
    mean_ratios=c()
    ratios = normalized_data[grepl("_ratio_",colnames(normalized_data))]
    for (ratio in ratios){
      mean = mean(ratio,na.rm=T)
      mean_ratios = c(mean_ratios,mean)
    }
    
    # Compute the normalized data
    i=1
    # mean_ref <- mean(as.numeric(as.vector(unlist(df[as.character(reference_sample)]))), na.rm = T)
    for (sample in colnames(df)[-1]){
      if(sample==reference_sample){
        normalized_data[paste0(reference_sample,"_norm")]<-reference_intensities
      }else{
        # normalized_data[paste0(sample,"_norm")]<-as.numeric(as.vector(unlist(df[as.character(sample)])))/normalized_data[paste0(sample,"_ratio_",reference_sample)]
        # normalized_data[paste0(sample,"_norm")]<-normalized_data[paste0(sample,"_ratio_",reference_sample)] * mean_ref
        normalized_data[paste0(sample,"_norm")]<-as.numeric(as.vector(unlist(df[as.character(sample)]))) * mean_ratios[i]
        i=i+1
      }
    }
  }
  
  return(normalized_data)
  
}

###################
#New Normalization#
###################

## 1) Normalization Factors ##

norm_factors <- function(subdf, reference_sample){
  # Initialize results dataframe
  ratio_names = paste0(colnames(subdf)[!colnames(subdf)%in%colnames(subdf)[1] & !colnames(subdf)%in%reference_sample],"_ratio_",reference_sample)
  norm_names = paste0(colnames(subdf)[-1],"_norm")
  normalized_data = data.frame(matrix(vector(),
                                      nrow(subdf),
                                      (length(ratio_names)+length(norm_names)+1),
                                      dimnames = list(c(), c("Id",ratio_names,norm_names) )),
                               check.names=FALSE)
  
  # Compute ratios
  normalized_data$Id = subdf[,1]
  reference_intensities = as.numeric(as.vector(unlist(subdf[as.character(reference_sample)])))
  for (sample in colnames(subdf)[-1]){
    if(sample!=reference_sample){
      normalized_data[paste0(sample,"_ratio_",reference_sample)]<-as.numeric(as.vector(unlist(subdf[as.character(sample)])))/reference_intensities
    }
  }
  
  # Compute median of ratios
  median_ratios=c()
  ratios = normalized_data[grepl("_ratio_",colnames(normalized_data))]
  for (ratio in ratios){
    median = median(ratio,na.rm=T)
    median_ratios = c(median_ratios,median)
  }
  
  return(median_ratios)
}

## 2) Create Dataframe ##

normdf <- function(raw_abundance, reference_sample){
  # Initialize results dataframe
  ratio_names = paste0(colnames(raw_abundance)[!colnames(raw_abundance)%in%colnames(raw_abundance)[1] & !colnames(raw_abundance)%in%reference_sample],"_ratio_",reference_sample)
  norm_names = paste0(colnames(raw_abundance)[-1],"_norm")
  normalized_df = data.frame(matrix(vector(),
                                    nrow(raw_abundance),
                                    (length(ratio_names)+length(norm_names)+1),
                                    dimnames = list(c(), c("Id",ratio_names,norm_names) )),
                             check.names=FALSE)
  
  # Compute ratios
  normalized_df$Id = raw_abundance[,1]
  reference_intensities = as.numeric(as.vector(unlist(raw_abundance[as.character(reference_sample)])))
  for (sample in colnames(raw_abundance)[-1]){
    if(sample!=reference_sample){
      normalized_df[paste0(sample,"_ratio_",reference_sample)]<-as.numeric(as.vector(unlist(raw_abundance[as.character(sample)])))/reference_intensities
    }
  }
  
  return(normalized_df)
}

## 3) Application ##

norm_app <- function(df, median_ratios, reference_intensities, normalized_data, raw_abundance){
  # Median removal
  i=1
  for (sample in colnames(df)[-1]){
    if(sample==reference_sample){
      normalized_data[paste0(reference_sample,"_norm")]<-raw_abundance[reference_sample]
    }else{
      normalized_data[paste0(sample,"_norm")]<-as.numeric(as.vector(unlist(df[as.character(sample)])))/median_ratios[i]
      i=i+1
    }
  }
  
  return(normalized_data)
}

## Aux Function ##

new_normalization <- function(sub_df, raw_abundance, reference_sample){
  med <- norm_factors(sub_df, reference_sample)
  
  normalized_df <- normdf(raw_abundance, reference_sample)
  
  reference_intensities = as.numeric(as.vector(unlist(sub_df[as.character(reference_sample)])))
  normalized_data <- norm_app(raw_abundance, med, reference_intensities, normalized_df, raw_abundance)
  
  return(normalized_data)
}


######
# CV #
######

# This function compute the median deviation for each sample of the quantification matrix.
compute_median_deviation <- function(intensity_matrix){
  
  # Input :
  # intensity_matrix : quantification matrix
  # Output :
  # matrix containing the median deviation for each intensity of each sample
  
  intensity_matrix$medians = rowMedians(data.matrix(intensity_matrix), na.rm=T)
  for(sample in colnames(intensity_matrix)[-ncol(intensity_matrix)]){
    intensity_matrix[sample]=(as.numeric(intensity_matrix[[sample]])-intensity_matrix$medians)/intensity_matrix$medians
  }
  
  return(intensity_matrix[,-which(names(intensity_matrix)=="medians")])
  
}

library(matrixStats)

# This function compute the median deviation for each sample of the quantification matrix. It requires the library "matrixStats".
compute_cv <- function(intensity_matrix){
  
  # Input :
  # intensity_matrix : quantification matrix
  # Output :
  # vector containing the CV for each protein
  
  medians = rowMedians(data.matrix(intensity_matrix), na.rm=T)
  sd = apply(intensity_matrix,1,sd)
  cv = sd/medians * 100
  
  return(cv)
  
}

##############
# Imputation #
##############
if(TRUE){
  # This function impute missing values with background noise using a gaussian model
  impute_background_noise_gaussian <- function(df){
    
    # Input :
    # df : quantification matrix
    # Output :
    # quantification matrix with imputed intensities
    
    # Remove "MCAR" annotations from the quantification matrix so numeric functions can work
    df_without_string = df
    df_without_string[df_without_string=="MCAR"] <- NA
    df_without_string = data.matrix(df_without_string)
    
    # Compute 1rst percentile (-> mean of the gaussian model)
    percentile = 0.01
    m = as.numeric(quantile(df_without_string[,-1],percentile, na.rm = T))
    
    # Compute sd between 1rst and 2nd percentile (-> sd of the gaussian model)
    min = m
    max = as.numeric(quantile(df_without_string[,-1],0.02, na.rm = T))
    intensities = as.numeric(unlist(df_without_string[,-1]))
    sd_exp = sd(intensities[intensities<=max & intensities>=min & !is.na(intensities)])
    
    # Create a gaussian distribution with the computed parameters
    n = sum(is.na(df))
    gaussian <- rnorm(2*n,m,sd_exp)
    
    # Pull missing values in gaussian distribution
    imputed_values = c()
    if(length(gaussian)<30){
      
      q1 = quantile(gaussian)[2]
      q3 = quantile(gaussian)[4]
      interval = gaussian[gaussian>q1 & gaussian<q3]
      imputed_values = c(sample(interval,n,replace=T))
      
    }else{
      
      gaussian = gaussian[gaussian>0]
      gaussian = sort(gaussian)
      
      ## Compute proportion of NA to pull in each interval
      interval_size = (max(gaussian)-min(gaussian))/21
      i=min(gaussian)
      while(i<max(gaussian)){
        interval = gaussian[gaussian>=i&gaussian<i+interval_size]
        percentage = length(interval)/length(gaussian)
        imputed_values = c(imputed_values,sample(interval,ceiling(percentage*n),replace=T))
        i = i+interval_size
      }
    }
    
    #Replace missing values in data
    df[is.na(df)] <- sample(imputed_values,length(df[is.na(df)]),replace=F)
    colnames(df)[2:ncol(df)] = paste0(colnames(df)[2:ncol(df)],"_imputed")
    
    results = list(m,sd_exp,df)
    
    return(results)
    
  }
  
  # This function impute missing values with background noise using a percentile.
  impute_background_noise_percentile <- function(df,per){
    
    # Input :
    # df : quantification matrix whom 1rst column contains protein id
    # per : chosen percentile
    # Output :
    # quantification matrix with imputed intensities
    
    for(col in colnames(df)[-1]){
      current_intensities = df[col]
      df_without_string = current_intensities
      df_without_string[df_without_string=="MCAR"] <- NA
      df_without_string = data.matrix(df_without_string)
      percentile = as.numeric(quantile(df_without_string,per, na.rm = T))
      current_intensities[is.na(current_intensities)] <- percentile
      df[col] = current_intensities
    }
    
    colnames(df)[2:ncol(df)] = paste0(colnames(df)[2:ncol(df)],"_imputed")
    return(df)
    
  }
  
  ########## MCAR Imputation ################
  
  # This function is dedicated to MCAR imputation.
  impute_mcar <- function(intensities,identification_type,model,threshold_mcar_obs,threshold_mcar_ms,knn_min_occurrences, data_f){
    
    # Input :
    # intensities : quantification matrix whom 1st column contains protein id
    # identification_type : matrix of identification type whom 1st column contains protein id
    # model : imputation model (KNN or none)
    # threshold_mcar_obs : minimum number of observations in the condition to be a MCAR
    # threshold_mcar_ms : minimum number of file identified by MS/MS in the condition to be a MCAR
    # knn_min_occurrences : minimum number of observations in the condition to be used as a nearest neighbor
    # filtered :filtered dataFrame for ranking
    # Output :
    # quantification matrix with imputed intensities
    
    colnames(intensities)[1] <- "Id"
    
    intensities_info = intensities
    identification_type_info = identification_type
    
    # Compute nb of NA by protein
    intensities_info$obs_sum = as.numeric(apply(intensities[-1],1,function(x) sum(!is.na(x))))
    intensities_info$na_sum = as.numeric(apply(intensities[-1],1,function(x) sum(is.na(x))))
    
    # Compute nb of MS/MS by protein
    identification_type_info$ms_sum = as.numeric(apply(identification_type_info[-1],1,function(x) sum(x=="By MS/MS")))
    
    # Order proteins by median intensity
    intensities_info$median = Biobase::rowMedians(as.matrix(intensities[-1]), na.rm=TRUE)
    intensities_info = intensities_info[order(intensities_info$median),]
    intensities_info$rank = seq(1:nrow(data_f))
    
    write.csv(intensities_info, "intensities_info.csv", row.names = F)
    write.csv(identification_type_info, "identification_type_info.csv", row.names = F)
    
    # Extract candidates for KNN
    putative_knn = subset(intensities_info,intensities_info$obs_sum>=knn_min_occurrences)
    putative_knn$cv = as.numeric(apply(putative_knn[,2:(ncol(putative_knn)-3)],1,function(x) sd(x,na.rm=T)/median(x,na.rm=T)))
    putative_knn$median = as.numeric(apply(putative_knn[,2:(ncol(putative_knn)-4)],1,function(x) median(x,na.rm=T)))
    
    # Find MCAR
    mcar = subset(identification_type_info,identification_type_info$ms_sum>=threshold_mcar_ms)
    mcar = subset(intensities_info,intensities_info$Id%in%mcar$Id)
    mcar = subset(mcar,(mcar$obs_sum>=threshold_mcar_obs & mcar$na_sum>0))
    
    write.csv(mcar, "mcar3.csv", row.names = F)
    
    # Impute MCAR
    if(nrow(mcar)>0){
      
      if(model=="none"){
        mcar[is.na(mcar)] <- "MCAR"
        
      }else{
        if(model=="knn"){
          
          # Number of nearest neighbours to select
          k = floor(sqrt(nrow(intensities)))
          
          for(row in 1:nrow(mcar)){
            
            current_mcar = mcar[row,]
            
            # Find knn
            current_rank = current_mcar$rank
            knn = putative_knn
            knn$distance = abs(putative_knn$rank - current_rank)
            knn = knn[order(knn$distance),]
            knn = knn[1:k,]
            
            # Informations about current MCAR
            aValues = as.numeric(current_mcar[,2:(ncol(current_mcar)-4)])
            nb_na_to_impute = sum(is.na(aValues)==T)
            median_protein_to_impute = median(aValues,na.rm=T)
            sd_knn = median(knn$cv,na.rm=T)*median_protein_to_impute # Expected sd after imputation
            q3 = as.numeric(quantile(aValues,na.rm=T)[4])
            d = q3 - median_protein_to_impute
            
            # Find intensities to frame sd_knn
            imputed = aValues
            rand =as.numeric(knn$median,na.rm=T)
            
            X = median_protein_to_impute
            imputed[is.na(imputed)] <- X
            
            while(sd(imputed)<sd_knn){
              imputed = aValues
              X = X+d
              imputed[is.na(imputed)] <- X
            }
            
            # Create the gaussian model
            if(X==median_protein_to_impute){
              gaussian = rnorm(10000, mean = median_protein_to_impute, sd = 0.05)
            }else{
              mean_X = (X + (X-d))/2
              gaussian = rnorm(10000, mean = mean_X, sd = 0.05)
            }
            
            # Replace NA by values sampled in the gaussian model
            rand = sample(gaussian,nb_na_to_impute,replace=T)
            imputed = aValues
            
            imputed[is.na(imputed)] <- rand
            
            mcar[row,2:ncol(intensities)] = imputed
          }
        }
      }
    }
    
    # Replace in the final table
    write("beg_end", "beg_end")
    write.csv(intensities, "test_intensities.csv", row.names = F)
    intensities = subset(intensities, !intensities[,1]%in%mcar[,1])
    intensities = rbind(intensities,mcar[,1:(ncol(mcar)-4)])
    write("beg_endd", "beg_endd")
    keep_col_names <- colnames(intensities)
    # intensities[is.nan(intensities)] <- NA
    intensities <- as.data.frame(lapply(intensities, function(x) ifelse(is.nan(x), NA, x)))
    write("beg_enddd", "beg_enddd")
    colnames(intensities) <- keep_col_names
    
    return(intensities)
  }
}


###########
# Heatmap #
###########

# This function create a heatmap with a hierarchical clustering of proteins
create_heatmap <- function(intensities,imputed,samples_names,group_ids,htitle){
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(ggdendro)
  library(patchwork)
  
  # Input:
  # intensities : quantification matrix where the first two columns are identifiers
  #               Example: Accession, Label, sample1, sample2, ...
  #
  # imputed : boolean matrix used as a missing-value mask
  #           Expected structure:
  #           Id, Label, sample1_mean, sample2_mean, ...
  #
  # samples_names : list of sample names without "_mean"
  #                 Example: K13D_1_C3_K13D, K13D_2_C3_K13D, ...
  #
  # group_ids : list of biological conditions
  #             Example: C3_K13D, C4_K13ND
  #
  # Output:
  # heatmap plot
  
  
  ############################
  # 1. Basic safety checks
  ############################
  
  # The heatmap matrix needs at least:
  # column 1 = Accession or protein identifier
  # column 2 = Label or gene name
  # column 3+ = intensity columns
  if (ncol(intensities) < 3) {
    stop("create_heatmap(): intensities must contain at least 3 columns: identifier, label, and intensity values.")
  }
  
  # The imputed mask should have the same number of protein rows as the intensity matrix.
  if (nrow(intensities) != nrow(imputed)) {
    stop("create_heatmap(): intensities and imputed do not have the same number of rows.")
  }
  
  # In your current setup:
  # intensities[[1]] is Accession
  # imputed[[1]] is named Id, but contains the same Accession values.
  #
  # This check makes sure the protein order is the same before any long-format reshaping.
  if (!identical(as.character(intensities[[1]]), as.character(imputed[[1]]))) {
    stop("create_heatmap(): protein order differs between intensities and imputed mask.")
  }
  
  
  ############################
  # 2. Build numeric matrix for clustering
  ############################
  
  # This replaces:
  # inputProteins = intensities[,3:ncol(intensities)]
  #
  # We keep only the sample intensity columns.
  inputProteins <- intensities[, 3:ncol(intensities), drop = FALSE]
  
  # Make sure all intensity columns are numeric.
  # This avoids hidden character columns causing strange behavior in dist/hclust/ggplot.
  inputProteins <- data.frame(
    lapply(inputProteins, as.numeric),
    check.names = FALSE
  )
  
  # Use gene labels as rownames for clustering.
  # make.names(..., unique = TRUE) prevents duplicate rowname problems.
  rownames(inputProteins) <- make.names(intensities[[2]], unique = TRUE)
  
  # Cluster proteins based on the intensity matrix.
  distMatrixProteins <- dist(inputProteins)
  hClustProteins <- hclust(distMatrixProteins)
  
  
  ############################
  # 3. Convert intensity matrix to long format
  ############################
  
  # These are the heatmap sample columns.
  # Example:
  # K13D_1_C3_K13D, K13D_2_C3_K13D, ...
  intensity_cols <- colnames(intensities)[3:ncol(intensities)]
  
  # This replaces the old manual construction:
  #
  # proteinRank = c()
  # for (i in 1:(ncol(inputProteins))){
  #   proteinRank = c(proteinRank,1:nrow(inputProteins))
  # }
  # proteinId = rep(intensities[,1],ncol(inputProteins))
  # gene_name = rep(make.names(intensities[,2],unique=T),ncol(inputProteins))
  #
  # The old version assumed a specific row order after pivoting.
  # This version keeps the protein identity attached to each row before pivoting.
  intensity_long <- intensities %>%
    dplyr::mutate(
      proteinRank = dplyr::row_number(),
      proteinId = .data[[colnames(intensities)[1]]],
      Gene_name = make.names(.data[[colnames(intensities)[2]]], unique = TRUE)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(intensity_cols),
      names_to = "Sample",
      values_to = "log_intensity"
    ) %>%
    dplyr::mutate(
      log_intensity = as.numeric(log_intensity)
    )
  
  
  ############################
  # 4. Convert imputed mask to long format
  ############################
  
  # Your imputed mask columns should contain "mean".
  # Example:
  # K13D_1_C3_K13D_mean, K13D_2_C3_K13D_mean, ...
  mask_cols <- grep("mean$", colnames(imputed), value = TRUE)
  
  if (length(mask_cols) == 0) {
    stop("create_heatmap(): imputed mask has no columns ending in '_mean'.")
  }
  
  # Check that the mask sample columns correspond to the heatmap sample columns.
  # The heatmap uses:
  # K13D_1_C3_K13D
  #
  # The mask uses:
  # K13D_1_C3_K13D_mean
  #
  # After removing "_mean", they should be identical.
  if (!identical(intensity_cols, sub("_mean$", "", mask_cols))) {
    stop(
      "create_heatmap(): heatmap sample columns and imputed mask columns do not match after removing '_mean'.\n",
      "Heatmap columns: ", paste(intensity_cols, collapse = ", "), "\n",
      "Mask columns: ", paste(mask_cols, collapse = ", ")
    )
  }
  
  # This replaces:
  #
  # pivot_long_imputed <- imputed %>%
  #   pivot_longer(cols=contains("mean"),
  #                names_to="condition",
  #                values_to="is_imputed")
  #
  # The important addition is proteinRank.
  # That gives us a stable way to join the mask to the intensity values.
  imputed_long <- imputed %>%
    dplyr::mutate(
      proteinRank = dplyr::row_number()
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(mask_cols),
      names_to = "Sample",
      values_to = "imputed"
    ) %>%
    dplyr::mutate(
      Sample = sub("_mean$", "", Sample),
      imputed = as.logical(imputed)
    )
  
  
  ############################
  # 5. Join intensity values and imputation mask
  ############################
  
  # This replaces:
  #
  # heatmapDF = data.frame(
  #   proteinId,
  #   gene_name,
  #   as.character(pivot_long_intensitiesDF$condition),
  #   pivot_long_intensitiesDF$intensity,
  #   proteinRank,
  #   pivot_long_imputed$is_imputed
  # )
  #
  # The old version put columns together by position.
  # The new version joins them by proteinRank and Sample.
  # This avoids misalignment when pivot_longer ordering is not what we assumed.
  heatmapDF <- intensity_long %>%
    dplyr::left_join(
      imputed_long[, c("proteinRank", "Sample", "imputed"), drop = FALSE],
      by = c("proteinRank", "Sample")
    )
  
  # If this happens, the mask did not match some heatmap cells.
  # That would make blank spaces unreliable, so we stop immediately.
  if (any(is.na(heatmapDF$imputed))) {
    stop("create_heatmap(): some heatmap cells do not have a matching imputation-mask value.")
  }
  
  # This is the line that creates blank spaces.
  # Imputed cells are removed from the plotted dataframe.
  # geom_raster() then leaves those positions blank.
  heatmapDF <- subset(heatmapDF, imputed == FALSE)
  
  if (nrow(heatmapDF) == 0) {
    stop("create_heatmap(): no non-imputed heatmap cells remain after filtering.")
  }
  
  
  ############################
  # 6. Order samples by group
  ############################
  
  # This replaces the original orderedSamples block.
  #
  # Important change:
  # I keep ordered_sample_levels separate from heatmapDF$orderedSamples.
  # In the old function, orderedSamples was first a sample list, then later reused
  # as a factor the length of the full heatmap dataframe.
  ordered_sample_levels <- c()
  
  group_ids <- sort(group_ids)
  
  for (group in group_ids) {
    ordered_sample_levels <- c(
      ordered_sample_levels,
      sort(grep(paste0("_", group), samples_names, value = TRUE))
    )
  }
  
  ordered_sample_levels <- unique(ordered_sample_levels)
  
  # Safety fallback:
  # If any intensity column was not captured by the group-based ordering,
  # append it at the end instead of creating NA x positions.
  missing_ordered_samples <- setdiff(intensity_cols, ordered_sample_levels)
  
  if (length(missing_ordered_samples) > 0) {
    ordered_sample_levels <- c(ordered_sample_levels, missing_ordered_samples)
  }
  
  heatmapDF$orderedSamples <- factor(
    heatmapDF$Sample,
    levels = ordered_sample_levels
  )
  
  heatmapDF$x <- as.numeric(heatmapDF$orderedSamples)
  
  if (any(is.na(heatmapDF$x))) {
    stop("create_heatmap(): some samples could not be assigned an x-axis position.")
  }
  
  
  ############################
  # 7. Order proteins according to clustering
  ############################
  
  # This keeps the original idea:
  # proteins are displayed in the order from hclust().
  heatmapDF$orderedProteins <- factor(
    heatmapDF$proteinRank,
    levels = hClustProteins$order
  )
  
  heatmapDF$y <- as.numeric(heatmapDF$orderedProteins)
  
  if (any(is.na(heatmapDF$y))) {
    stop("create_heatmap(): some proteins could not be assigned a y-axis position.")
  }
  
  # Labels in clustered order.
  # This is safer than using unique(heatmapDF$Gene_name), because heatmapDF has
  # multiple rows per protein and some cells have been removed.
  y_axis_labels <- make.names(intensities[[2]], unique = TRUE)[hClustProteins$order]
  
  
  ############################
  # 8. Adjust color limits
  ############################
  
  # These were in the original function.
  # They are currently not used directly in scale_fill_gradientn(),
  # but I am keeping them in case you use them later.
  minRatio <- min(heatmapDF$log_intensity, na.rm = TRUE)
  maxRatio <- max(heatmapDF$log_intensity, na.rm = TRUE)
  maxabsRatio <- max(abs(heatmapDF$log_intensity), na.rm = TRUE)
  
  
  ############################
  # 9. Make heatmap plot
  ############################
  
  hcPlot <- ggplot(data = heatmapDF) +
    geom_raster(mapping = aes(x = x, y = y, fill = log_intensity)) +
    scale_fill_gradientn(colours = c("green", "yellow", "red")) +
    ggtitle("") +
    scale_x_continuous(
      name = "",
      breaks = seq_along(ordered_sample_levels),
      labels = ordered_sample_levels,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      name = "",
      breaks = seq_along(hClustProteins$order),
      labels = y_axis_labels,
      expand = c(0, 0),
      position = "right"
    ) +
    labs(title = htitle) +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 30),
      axis.title = element_text(size = 25),
      axis.text.x = element_text(size = 24, angle = 90),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)
    )
  
  
  ############################
  # 10. Make protein dendrogram
  ############################
  
  dendrogramProteins <- dendro_data(hClustProteins, type = "rectangle")
  dendrogramDataProteins <- segment(dendrogramProteins)
  
  dendrogramProteinsPlot <- ggplot() +
    geom_segment(
      data = dendrogramDataProteins,
      mapping = aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    coord_flip(clip = "off") +
    scale_y_reverse(expand = c(0.5, 0)) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  
  
  ############################
  # 11. Combine dendrogram and heatmap
  ############################
  
  result <- dendrogramProteinsPlot +
    plot_spacer() +
    hcPlot +
    plot_layout(
      widths = c(4, -1.1, 4.5),
      guides = "collect"
    )
  
  # These were in your original function.
  result$heights[7] <- unit(0.2, "null")
  result$widths[5] <- unit(0.2, "null")
  
  return(result)
}

########################
#Heatmap w/o clustering#
########################

# This function create a heatmap with a hierarchical clustering of proteins
create_heatmap_2 <- function(intensities,imputed,samples_names,group_ids){
  
  # Input :
  # intensities : quantification matrix whom 1rst column contains protein id
  # imputed : boolean matrix whom 1rst column contains protein id
  # samples_names : list of sample names
  # group_ids : list of biological conditions
  # Output :
  # heatmap plot
  
  inputProteins = intensities[,3:ncol(intensities)]
  rownames(inputProteins) = make.names(intensities$Label,unique=T)
  distMatrixProteins <- dist(inputProteins)
  hClustProteins <- hclust(distMatrixProteins)
  
  # Create data frame
  proteinRank = c()
  for (i in 1:(ncol(inputProteins))){
    proteinRank = c(proteinRank,1:nrow(inputProteins))
  }
  proteinId = rep(intensities[,1],ncol(inputProteins))
  gene_name = rep(make.names(intensities[,2],unique=T),ncol(inputProteins))
  intensitiesDF = melt(inputProteins)
  imputed = melt(imputed,id.vars = 1)
  heatmapDF = cbind(proteinId,gene_name,intensitiesDF,proteinRank,imputed$value)
  colnames(heatmapDF) = c("Id","Gene_name","Sample","log_intensity","proteinRank","imputed")
  
  heatmapDF=subset(heatmapDF,heatmapDF$imputed==F)
  
  # Order sample by group
  orderedSamples = c()
  group_ids = sort(group_ids)
  for (group in group_ids){
    orderedSamples = c(orderedSamples,sort(grep(paste0("_",group),samples_names,value=T)))
  }
  heatmapDF$orderedSamples <- factor(heatmapDF$Sample, levels = unique(orderedSamples))
  heatmapDF$x <- as.numeric(heatmapDF$orderedSamples)
  
  # Order proteins according to clustering results
  heatmapDF$orderedProteins <- factor(heatmapDF$proteinRank, levels = hClustProteins$order)
  heatmapDF$y <- as.numeric(heatmapDF$orderedProteins)
  
  # Adjust colors
  minRatio <- min(heatmapDF$log_intensity)
  maxRatio <- max(heatmapDF$log_intensity)
  maxabsRatio <- max(abs(heatmapDF$log_intensity))
  
  # Transforming Heatmap DataFrame
  # heatmapDF <- t(heatmapDF)
  # write.csv(heatmapDF, "heatmapDF")
  
  # Make plot
  hcPlot <- ggplot(data = heatmapDF) +
    geom_raster(mapping = aes(x = y, y = x, fill = log_intensity)) +
    scale_fill_gradientn(colours=c("green","yellow","red")) +
    #scale_fill_scico(palette = "berlin") +
    ggtitle("")+
    scale_y_continuous(name = "", breaks = 1:length(orderedSamples), labels = orderedSamples, expand = c(0, 0), position = "right") +
    scale_x_continuous(name = "", breaks = unique(heatmapDF$y), labels = unique(heatmapDF$Gene_name), expand = c(0, 0)) +
    theme(#axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title=element_text(size=30),
      axis.title=element_text(size=25),
      axis.text.y=element_text(size=24),
      axis.text.x=element_text(angle = 90),
      legend.text=element_text(size=16),
      legend.title=element_text(size=18))
  
  return(hcPlot)
  
  # No dendrogram needed
  
  # Make grobs from plots
  # matrixGrob <- ggplotGrob(hcPlot)
  
  # Assemble
  # bottom
}

#######
# ROC #
#######

# This function compute the sensibility and the specificity varying ratio and pvalue thresholds
compute_ROC <- function(statistic_table,n,ratioVar,pvalVar,variants){
  
  # Input :
  # statistic_table : table containing results of statistical analysis and whom 1rst column corresponds to protein ids
  # n : number of tests
  # ratioVar : ratio thresholds to test
  # pvalVar : pval thresholds to test
  # variants : list of ids of expected variants
  # Output :
  # heatmap plot
  
  mat <- matrix(ncol = 4, nrow = (n*n))
  dimnames(mat)[[2]] <- c("FDP", "Sensitivity", "ratio_threshold", "pval_threshold")
  l <- 0
  for (i in 1:n) {
    Z <- ratioVar[i]
    for (j in 1:n) {
      p <- pvalVar[j]
      
      FP <- nrow(statistic_table[(abs(statistic_table$ratio) >= Z) & (statistic_table$pval <= p) & !(statistic_data$Accession %in% variants),])
      TP <- nrow(statistic_table[(abs(statistic_table$ratio) >= Z) & (statistic_table$pval <= p) & (statistic_data$Accession %in% variants),])
      TN <- nrow(statistic_table[((abs(statistic_table$ratio) < Z) | (statistic_table$pval > p)) & !(statistic_data$Accession %in% variants),])
      FN <- nrow(statistic_table[((abs(statistic_table$ratio) < Z) | (statistic_table$pval > p)) & (statistic_data$Accession %in% variants),])
      
      speci <- (FP/(FP+TP))*100
      sensi <- (TP/(TP+FN))*100
      mat[j+l,] <- c(speci, sensi, Z, p)
    }
    l <- n+l
  }
  gtab <- as.data.frame(mat)
  return(gtab)
  
}