
#################
# Normalization #
#################

normalize <- function(df, reference_sample){

  # Initialize dataframe
  ratio_names = paste0(colnames(df)[!colnames(df)%in%colnames(df)[1] & !colnames(df)%in%reference_sample],"_ratio_",reference_sample)
  norm_names = paste0(colnames(df)[-1],"_norm")
  normalized_data = data.frame(matrix(vector(),
                                      nrow(df),
                                      (length(ratio_names)+length(norm_names)+1),
                                      dimnames = list(c(), c("Id",ratio_names,norm_names) )),
                               check.names=FALSE)

  # Compute ratios
  normalized_data$Id = df[,1]
  reference_intensities = as.numeric(as.vector(unlist(df[as.character(reference_sample)])))
  for (sample in colnames(df)[-1]){
    if(sample!=reference_sample){
      normalized_data[paste0(sample,"_ratio_",reference_sample)]<-as.numeric(as.vector(unlist(df[as.character(sample)])))/reference_intensities
    }
  }

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

  return(normalized_data)

}


######
# CV #
######

compute_median_deviation <- function(intensity_matrix){

  intensity_matrix$medians = rowMedians(data.matrix(intensity_matrix), na.rm=T)
  for(sample in colnames(intensity_matrix)[-ncol(intensity_matrix)]){
    intensity_matrix[sample]=(as.numeric(intensity_matrix[[sample]])-intensity_matrix$medians)/intensity_matrix$medians
  }

  return(intensity_matrix[,-which(names(intensity_matrix)=="medians")])

}

library(matrixStats)

compute_cv <- function(intensity_matrix){

  #means = rowMeans(as.matrix(intensity_matrix), na.rm=T)
  medians = rowMedians(data.matrix(intensity_matrix), na.rm=T)
  sd = apply(intensity_matrix,1,sd)
  cv = sd/medians * 100

  return(cv)

}

##############
# Imputation #
##############

impute_background_noise_gaussian <- function(df){

  # Remove MCAR annotations so numeric functions can work
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


impute_background_noise_percentile <- function(df,per){

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

impute_mcar <- function(intensities,identification_type,model,threshold_mcar_obs,threshold_mcar_ms,knn_min_occurrences){

  intensities_info = intensities
  identification_type_info = identification_type

  # Compute nb of NA
  intensities_info$obs_sum = as.numeric(apply(intensities[-1],1,function(x) sum(!is.na(x))))
  intensities_info$na_sum = as.numeric(apply(intensities[-1],1,function(x) sum(is.na(x))))

  # Compute nb of MS/MS
  identification_type_info$ms_sum = as.numeric(apply(identification_type_info[-1],1,function(x) sum(x=="By MS/MS")))

  # Order proteins by median intensity
  intensities_info$median = Biobase::rowMedians(as.matrix(intensities[-1]), na.rm=TRUE)
  intensities_info = intensities_info[order(intensities_info$median),]
  intensities_info$rank = seq(1:nrow(filtered_data))

  # Extract candidates for KNN
  putative_knn = subset(intensities_info,intensities_info$obs_sum>=knn_min_occurrences)
  putative_knn$cv = as.numeric(apply(putative_knn[,2:(ncol(putative_knn)-3)],1,function(x) sd(x,na.rm=T)/median(x,na.rm=T)))
  putative_knn$median = as.numeric(apply(putative_knn[,2:(ncol(putative_knn)-4)],1,function(x) median(x,na.rm=T)))

  # Find MCAR
  mcar = subset(identification_type_info,identification_type_info$ms_sum>=threshold_mcar_ms)
  mcar = subset(intensities_info,intensities_info$`filtered_data$Id`%in%mcar$Id)
  mcar = subset(mcar,(mcar$obs_sum>=threshold_mcar_obs & mcar$na_sum>0))

  if(nrow(mcar)>0){

    if(model=="none"){
      mcar[is.na(mcar)] <- "MCAR"

    }else{
      if(model=="knn"){

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

          # X = median_protein_to_impute
          # imputed[is.na(imputed)] <- X
          #
          # while(sd(imputed)<sd_knn){
          #   imputed = aValues
          #   X = X+d
          #   imputed[is.na(imputed)] <- X
          # }
          #
          # # Create the gaussian model
          # if(X==median_protein_to_impute){
          #   gaussian = rnorm(10000, mean = median_protein_to_impute, sd = 0.05)
          # }else{
          #   mean_X = (X + (X-d))/2
          #   gaussian = rnorm(10000, mean = mean_X, sd = 0.05)
          # }
          #
          # # Impute NA
          # rand = sample(gaussian,nb_na_to_impute,replace=T)
          # imputed = aValues


          imputed[is.na(imputed)] <- rand

          mcar[row,2:ncol(intensities)] = imputed

        }

      }
    }

  }

  # Replace in the final table
  intensities = subset(intensities, !intensities[,1]%in%mcar[,1])
  intensities = rbind(intensities,mcar[,1:(ncol(mcar)-4)])
  intensities[intensities==NaN] <- NA

  return(intensities)
}


###########
# Heatmap #
###########

create_heatmap <- function(intensities,imputed,samples_names,group_ids){

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
  heatmapDF$orderedSamples <- factor(heatmapDF$Sample, levels = orderedSamples)
  heatmapDF$x <- as.numeric(heatmapDF$orderedSamples)

  # Order proteins according to clustering results
  heatmapDF$orderedProteins <- factor(heatmapDF$proteinRank, levels = hClustProteins$order)
  heatmapDF$y <- as.numeric(heatmapDF$orderedProteins)

  # Adjust colors
  minRatio <- min(heatmapDF$log_intensity)
  maxRatio <- max(heatmapDF$log_intensity)
  maxabsRatio <- max(abs(heatmapDF$log_intensity))

  # Make plot

  hcPlot <- ggplot(data = heatmapDF) +
    geom_raster(mapping = aes(x = x, y = y, fill = log_intensity)) +
    scale_fill_gradientn(colours=c("green","yellow","red")) +
    #scale_fill_scico(palette = "berlin") +
    ggtitle("")+
    scale_x_continuous(name = "", breaks = 1:length(orderedSamples), labels = orderedSamples, expand = c(0, 0)) +
    scale_y_continuous(name = "", breaks = unique(heatmapDF$y), labels = unique(heatmapDF$Gene_name), expand = c(0, 0), position = "right") +
    theme(#axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title=element_text(size=30),
          axis.title=element_text(size=25),
          axis.text.x=element_text(angle=90, size=24),
          #axis.text.y=element_text(size=12),
          legend.text=element_text(size=16),
          legend.title=element_text(size=18))

  dendrogramProteins <- dendro_data(hClustProteins, type="rectangle")
  dendrogramDataProteins <- segment(dendrogramProteins)
  dendrogramProteinsPlot <- ggplot() +
    geom_segment(data = dendrogramDataProteins, mapping = aes(x=x, y=y, xend=xend, yend=yend)) +
    #geom_text(data=label(dendrogramProteins), aes(x=x, y=y, label=label, hjust=0), size=3) +
    coord_flip(clip = 'off') +
    scale_y_reverse(expand = c(0.5, 0)) +
    #expand_limits(x = c(10, 10)) +
    #scale_x_continuous(expand = c(0, 0), limits = c(0.5, nrow(intensities) + 0.5)) +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank())

  # Make grobs from plots
  matrixGrob <- ggplotGrob(hcPlot)
  dendroProteinGrob <- ggplotGrob(dendrogramProteinsPlot)

  # Assemble
  bottomGrob <- cbind(dendroProteinGrob[, 5], matrixGrob[, 5:ncol(matrixGrob)], size = "last")
  bottomGrob <- cbind(matrixGrob[, 1:4], bottomGrob, size = "last")

  result <- rbind(bottomGrob[7:nrow(bottomGrob), ], size = "last")
  result <- rbind(bottomGrob[1:6, ], result, size = "last")

  # Adjust sizes
  result$heights[7] <- unit(0.2, "null")
  result$widths[5] <- unit(0.2, "null")

  return(result)
}

#######
# ROC #
#######

compute_ROC <- function(statistic_table,n,ratioVar,pvalVar,variants){

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
