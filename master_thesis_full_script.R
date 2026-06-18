# Master thesis script
# Author: Nele Noack

---------------------------##### Field Data Analysis ##### ---------------------------------------------------------------------
library(data.table)
# Load data
input <- if (file.exists("fielddata.csv")) {
  "fielddata.csv"
} else {
  ".../fielddata.csv"
}
fielddata <- fread(input)

# extract the number of trees for each plot
number_trees <- fielddata[,.(nr_trees = .N), by = plot]

# determine species richness, abundance, Shannon and Simpson over all plots
# for each plot
total_species_richness <- fielddata[,.(species_richness = uniqueN(species_name))]

overall_species <- fielddata[,.(species_count = .N), by = .(species_name)]
overall_abundance <- fielddata[,.(species_count = sum(species_count)), by = .(species_name)]
overall_shannon <- fielddata[, .(Shannon = sum(-(species_count/sum(species_count))*log(species_count/sum(species_count))))]

# Determine species diversity indices per plot
# determine number of trees per species
nr_trees_species <- fielddata[,.(species_count = .N), by = .(plot, species_name)]

# Species richness
species_richness <- fielddata[,.(species_richness = uniqueN(species_name)), by = plot]
species_richness

# Shannon index: sum(nr_trees_species/total_nr_trees * log (nr_trees_species/total_nr_trees))
shannon_index <- nr_trees_species[, .(Shannon = sum(-(species_count/sum(species_count))*log(species_count/sum(species_count)))) , by = plot]

# Simpson index: sum((nr_trees_species/total_nr_trees)²)
simpson_index <- nr_trees_species[, .(Simpson = sum((species_count/sum(species_count))**2)), by = plot]
inverse_simpson <- nr_trees_species[, .(Inv_Simpson = 1/sum((species_count/sum(species_count))**2)), by = plot]

# merge the indices into one table
merged_diversity <- merge(species_richness, shannon_index, by="plot")
merged_diversity <- merge(merged_diversity, inverse_simpson, by="plot")
merged_diversity <- merge(merged_diversity, number_trees, by="plot")

# load coordinate data and merge with diversity data
input_c <- if (file.exists("Shagayu_Field_Data_coordinates_reviewed.csv")) {
  "Shagayu_Field_Data_coordinates_reviewed.csv"
} else {
  "../coordinates.csv"
}
coordinates <- fread(input_c)
diversity_coordinates <- merge(coordinates, merged_diversity, by="plot")
write.csv(diversity_coordinates, file = 'species_diversity.csv', row.names = FALSE)


------------------------##### Satellite image preprocessing #####------------------------------------------------------------------
library(terra)
library(sf)

# load sentinel-1 and sentinel-2 raster images
sen1 <- rast("Sentinel1_reprojectedArc1960.tif")
sen2 <- rast("Sentinel2_reprojectedArc1960.tif")
# stack sentinel-2 and sentinel-1
stacked_sen <- c(sen2, sen1)

### calculate the vegetation indices
# NDVI: with Sentinel-2 B8 (NIR - here index 3) and B4 (Red here index 7)
ndvi <- ((sen2[[7]] - sen2[[3]])/(sen2[[7]] + sen2[[3]]))
# NDWI: sentinel-2 B8 (NIR) and B11 (SWIR1)
ndwi <- ((sen2[[7]] - sen2[[9]])/(sen2[[7]] + sen2[[9]]))
# NDPI	Normalized difference phenology index	(NIR - (0.74 × Red + 0.26 × SWIR1)) / (NIR + (0.74 × Red + 0.26 × SWIR1))
ndpi <- ((sen2[[7]] - (0.74*sen2[[3]] + 0.26*sen2[[9]]))/(sen2[[7]] + (0.74*sen2[[3]] + 0.26*sen2[[9]])))
# EVI	Enhanced vegetation index	2.5 × (NIR - Red) / (NIR + 6 × Red - 7.5 × Blue + 1)
evi <- (2.5*((sen2[[7]] - sen2[[3]])/(sen2[[7]] + 6*sen2[[3]] - (7.5*sen2[[1]]) + 1)))
# EVI2	2-band Enhanced Vegetation Index	2.5x(NIR-Red) / (NIR + 2.4Red +1)
evi2 <- (2.5*((sen2[[7]] - sen2[[3]])/(sen2[[7]] + 2.4*sen2[[3]] + 1)))
# NBR	Normalized burn ratio index	(NIR - SWIR2) / (NIR + SWIR2)
nbr <- ((sen2[[7]] - sen2[[10]])/((sen2[[7]] + sen2[[10]])))
# NBR2	Normalized burn ratio index 2	(SWIR1 - SWIR2) / (SWIR1 + SWIR2)
nbr2 <- ((sen2[[9]] - sen2[[10]])/((sen2[[9]] + sen2[[10]])))
# NIRv	Near-infrared reflectance of terrestrial vegetation	NIR × NDVI
nirv <- (sen2[[7]] * ndvi)
# PI Phenology Index	NDVI²-NDII² with NDII = (NIR - SWIR1)/(NIR + SWIR)
phen_ind <- (ndvi**2 - ndmi**2)
# CLre Red-Edge Chlorophyll index	(RedEdge3/RedEdge1) - 1
clre <- ((sen2[[6]]/sen2[[4]]) - 1)
# CCI	Chlorophyll Carotenoid index	(Green - Red)/(Green + Red)
cci <- (sen2[[2]] - sen2[[3]])/(sen2[[2]] + sen2[[3]])
# IRECI	Inverted Red-Edge Chlorophyll Index	(NIR-Red)/(RE1/RE2)
ireci <- (sen2[[7]] - sen2[[3]])/(sen2[[4]]/sen2[[5]])
# NDI45	NDI for band 4 and 5	(RE1 - Red)/(RE1 + Red)
ndi45 <- (sen2[[4]] - sen2[[3]])/(sen2[[4]] + sen2[[3]])
# SR Simple Ratio Index	(NIR/Red)
sr <- (sen2[[7]]/sen2[[3]])
# SAVI	Soil Adjusted Vegetation Index	1.5x(NIR1 - Red)/(NIR1 + Red + 0.5)
savi <- (1.5 * (sen2[[7]] - sen2[[3]])/(sen2[[7]] + sen2[[3]] + 0.5))

# stack the layers in one raster 
stacked_sentinel_complete <- c(stacked_sen, ndvi, ndmi, ndpi, evi, evi2, nbr, 
                               nbr2, nirv, phen_ind, clre, cci, ireci, ndi45, sr, savi)
# rename the raster layers
names(stacked_sentinel_complete) <- c("Sentinel2_B2", "Sentinel2_B3", "Sentinel2_B4","Sentinel2_B5", "Sentinel2_B6",
                            "Sentinel2_B7", "Sentinel2_B8", "Sentinel2_B8A", "Sentinel2_B11", "Sentinel2_B12",
                            "Sentinel1_VV", "Sentinel1_VH", "NDVI", "NDWI", "NDPI", "EVI", "EVI2", "NBR", "NBR2",
                            "NIRv", "PI", "Clre", "CCI", "IRECI", "NDI45", "SR", "SAVI")
# clip raster using clipping extent of Shagayu Forest reserve
clipped_img <- rast("clip.tif")
cropped_stack <- crop(stacked_sentinel_complete, clipped_img)
# mask out non-forest pixels from the landcover mask
veg_mask <- rast("treecover_mask_2024.tif")
masked_sentinel_raster <- mask(x=cropped_stack, mask=veg_mask)

writeRaster(masked_sentinel_raster, "stacked_sentinel_complete.tif")

--------------------------------##### Spectral Heterogeneity Metrics #####------------------------------------------------------
# load the sentinel stack
sentinel_raster <- rast("stacked_sentinel27_complete.tif")


# normalize the data - use stretch to 0 to 255
sentinel_raster_8bit <- terra::stretch(sentinel_raster, minv=0, maxv=255)

# determine mean, standard deviation and coefficient of variation using the terra package
# moving window with focal() function - w is the window size, fun is the function, sd is standard deviation
# Same for mean value
mean_rast <- focal(sentinel_raster_8bit, w=3, fun="mean")
writeRaster(mean_all_test, "mean_allLayers_3.tif")

stdv_rast <- focal(sentinel_raster_8bit, w=3, fun = "sd")
writeRaster(stdv_all_test, "sdtv_allLayers_3.tif")

# Coefficient of variation as ratio between stdv and mean
cv_rast <- stdv_rast / mean_rast
writeRaster(cv_rast, "cv_allLayers_3.tif")


# Rao's Q with the rasterdiv package
library(rasterdiv)
library(rasterVis)
library(RColorBrewer)

# rasterdiv uses SpatRaster data format -> use rast()
sentinel_raster <- rast("stacked_sentinel27_complete.tif")

# stretch to 8-bit raster for rasterdiv applications and round for discrete data
sentinel_raster_8bit <- round(terra::stretch(masked_sentinel_raster, minv=0, maxv=255))

# run parametric Rao for each band/layer using the distance weight (alpha)
# here this was run for alpha = 1, alpha = 2 and alpha = 5 separately
# window is the window size, na.tolerance determines what proportion of pixel values can be empty
parao_lyr <- paRao(sentinel_raster_8bit[[1]], window = 3, alpha= 1, na.tolerance = 0.5)
paRao_layered_raster <- rast(parao_lyr_1[[1]])

for (lyr in 2:27){
  paRao_per_layer <- paRao(sentinel_raster_8bit[[lyr]], window = 3, alpha= 1, na.tolerance = 0.5)
  paRao_layered_raster <- c(paRao_layered_raster, rast(paRao_per_layer[[1]]))
  
}
terra::writeRaster(paRao_layered_raster, "paRao_allLayers.tif")


### Spectral Species using the biodivMapR package
library(preprocS2)
library(biodivMapR)
input_file_path_sent <- file.path("biodivMapR", "stacked_sentinel27.tif")
sentinel_raster <- rast("stacked_sentinel27.tif")
mask_path <- file.path("biodivMapR", "treecover_mask_2024.tif")
treemask <- rast("treecover_mask_2024.tif")

# 1- define biodivMapR output directory 
output_dir <- './biodivMapR/Results'
dir.create(path = output_dir, showWarnings = F, recursive = T)
output_dir_biodivMapR_test <- file.path(output_dir, 'fullbands_50')
dir.create(output_dir_biodivMapR_test, showWarnings = F, recursive = T)

# 2- define parameters for biodivMapR
window_size <- 3     # window size for computation of spectral diversity
nb_clusters <- 50     # number of clusters (spectral species) - here 20, 50 and 100 were used

# 3- define path for intermediate variables to be saved
# - information related to kmeans clustering defining spectral species
Kmeans_info_save <- file.path(output_dir_biodivMapR_test,'Kmeans_info.RData')
# - information related to beta diversity mapping (BC dissimilarity + PCoA)
Beta_info_save <- file.path(output_dir_biodivMapR_test,'Beta_info.RData')

# 4- adjust parameters for multithread & computational efficiency
nbCPU <- 3           # nb of threads for parallel processing

# 5- apply biodivMapR - k-means clustering and calculation of the spectral species diversity indices for 3x3 window
selected_bands <- c(1:27)
options <- set_options_biodivMapR(fun = 'biodivMapR_full')
options$nb_clusters <- nb_clusters
options$maxRows <- 1000
options$alpha_metrics <- c('richness', 'shannon','simpson')
ab_info_SPCA <- biodivMapR_full(input_raster_path = input_file_path_sent, 
                                input_mask_path = mask_path,
                                output_dir = output_dir_biodivMapR_test, 
                                selected_bands = selected_bands, 
                                window_size = window_size, 
                                Kmeans_info_save = Kmeans_info_save,
                                Beta_info_save = Beta_info_save, 
                                nbCPU = nbCPU, options = options)

------------------------------------##### Correlation Analysis #####------------------------------------------------------------
library("PerformanceAnalytics")
library(corrplot)
library("Hmisc")
library("dplyr")
library(ggplot2)    # Graphics library
library(sf)         # Spatial data types and handling
library(mapview)    # Visualize spatial data
library(spdep)      # Diagnosing spatial dependence
library(spatialreg) # Spatial lag and spatial error model
library(gstat)
library(sp)
library(car)

#### correlation analysis for the field data
# load the species diversity data
input <- if (file.exists("fielddata_diversity_reviewed.csv")) {
  "fielddata_diversity_reviewed.csv"
} else {
  ".../fielddata_diversity_reviewed.csv"
}
fielddata_indices <- read.csv(input)

#  split into training and test data
training_fielddata  <- fielddata_indices %>% dplyr::sample_frac(0.7)
write.csv(training_fielddata, file = 'training_fielddata_1.csv', row.names = FALSE)
training_fielddata$ID <- seq.int(nrow(training_fielddata))

testing_fielddata   <- dplyr::anti_join(fielddata_indices, training_fielddata, by = 'Plot')
write.csv(testing_fielddata, file = 'testing_fielddata_1.csv', row.names = FALSE)
testing_fielddata$ID <- seq.int(nrow(testing_fielddata))

# create spat Vector data from training data
training_pts_vect <- vect(training_fielddata, crs="EPSG:21037", geom=c("X", "Y"))


# spatial autocorrelation test
##### ---- Semivariogram ----- #####
test_training_autoc <- training_fielddata 
coordinates(test_training_autoc) <-c ("X","Y")
class(test_training_autoc)
hscat(test_training_autoc$Shannon ~ 1, test_training_autoc, (0:15) * 300) 

shannon_var <- variogram(test_training_autoc$Shannon ~ 1, test_training_autoc)
plot(shannon_var, col = "black", 
     main = "Variogram for Shannon's H comparing points from training data", xlab = "distance in meters")

hscat(test_training_autoc$species_richness ~ 1, test_training_autoc, (0:15) * 150) 
richness_var <- variogram(test_training_autoc$species_richness ~ 1, test_training_autoc, width=500)
plot(richness_var, col = "black", 
     main = "Variogram for species richness comparing points from training data", xlab = "distance in meters")

hscat(test_training_autoc$Inv_Simpson ~ 1, test_training_autoc, (0:15) * 400) 
simpson_var <- variogram(test_training_autoc$Inv_Simpson ~ 1, test_training_autoc, width= 400)
plot(simpson_var, col = "black", 
     main = "Variogram for Inverse Simpson Index comparing points from training data", xlab = "distance in meters")


#### ---- correlation between spectral metrics and species diversity metrics ---- ######
matrix_fielddata <- data.matrix(training_fielddata[,4:7])
chart.Correlation(matrix_fielddata, histogram=TRUE, pch=19, method="spearman")


### ----- Correlation with spectral heterogeneity metrics ----- ####
### here examplified using the CV but same approach for mean, standard deviation and Rao's Q
cv_rast <- rast("cv_allLayers.tif")

# extract values at training points from raster data
cv_points_extract <- extract(cv_rast, training_pts_vect, raw = T, ID = T)
cv_points_df <- data.frame(cv_points_extract)
# remove NoData values
cv_points_df <- na.omit(cv_points_df)

# test for normality of selected bands
qqnorm(cv_points_df$Sentinel2_B2, main = "a) Q-Q Plot of CV (Sentinel-2 Band 2)")
qqline(cv_points_df$Sentinel2_B2, col = "steelblue", lwd = 2)
shapiro.test(cv_points_df$Sentinel2_B2)

# correlation between (CV) bands
correlation_cv <- rcorr(cv_points_extract, type = "spearman")
correlation_cv$P  # p-value

# merge field species diversity training data and spectral heterogeneity training data into one data frame / matrix
df_cv_fielddata <- merge(cv_points_df, training_fielddata, by="ID")
# runc spearman rank correlation 
cv_corr_species <- rcorr(matrix_fielddata, cv_points_extract, type = "spearman")
# plot correlation
corrplot(cv_corr_species$r, method="number", type="lower", order="hclust", 
         p.mat = cv_corr_species$P, sig.level = 0.05, insig = "blank", title = "CV - Secies diversity Correlogram from rcorr function - w/ significance level")
# save the correlation coefficient (r) and the significanc score (p)
write.csv(cv_corr_species$r, file = 'correlation_cv_diversity.csv', row.names = FALSE)
write.csv(cv_corr_species$P, file = 'correlationP_cv_diversity.csv', row.names = FALSE)

# plot cv(CCI) and species richness
sp_cv <- cor.test(df_cv_fielddata$CCI, df_cv_fielddata$species_richness, method = "spearman")
rho_cv <- round(sp_cv$estimate, 3)
pval_cv <- signif(sp_cv$p.value, 3)

ggplot(df_test_cv, aes(CCI, species_richness)) +
  geom_point(color = "steelblue") +
  geom_smooth(method = "loess", color = "red") +
  labs(title = "a) Coefficient of Variation (CCI) ~ Richness") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman ρ = ", rho_cv, "\n",
                          "p = ", pval_cv),
           hjust = 1.1, vjust = 1.5,
           size = 5)
# the same was repeated for mean, standard deviation, rao's q (with distance weight 1, 2, 5)


#### correlation with Spectral Species Diversity ####
# load spectral species diversity rasters
spectral_richness_rast <- rast("richness_mean_20.tiff")
spectral_shannon_rast <- rast("shannon_mean_20.tiff")
spectral_simpson_rast <- rast("simpson_mean_20.tiff")
# stack spectral species diversity rasters into one
spectral_species_div <- c(spectral_richness_rast, spectral_shannon_rast, spectral_simpson_rast)
# extract values from training points
ssdiv_points_extract <- extract(spectral_species_div50, training_pts_vect, raw = T, ID = F)

ssdiv_points_df <- data.frame(ssdiv_points_extract)
# q-q plot
qqnorm(ssdiv_points_df$mean.shannon, main = "f) Q-Q Plot of Spectral Species Shannon")
qqline(ssdiv_points_df$mean.shannon, col = "steelblue", lwd = 2)
# normality test
shapiro.test(ssdiv_points_df$mean.shannon)
# correlation analysis of spectral species diversity metrics and tree species diversity indices
ssdiv_corr_species <- rcorr(matrix_fielddata, ssdiv_points_extract, type = "spearman")
ssdiv_corr_species$r
corrplot(ssdiv_corr_species$r, method="number", type="lower", order="hclust", 
         p.mat = ssdiv_corr_species$P, sig.level = 0.05, insig = "blank", title = "RaoQ 2 - Secies diversity Correlogram from rcorr function - w/ significance level")

chart.Correlation(data.matrix(ssdiv_points_extract), method = "spearman", histogram=TRUE, pch=19)
write.csv(ssdiv_corr_species$r, file = 'correlation_ssdiv20.csv', row.names = FALSE)
write.csv(ssdiv_corr_species$P, file = 'correlationP_ssdiv20.csv', row.names = FALSE)


---------------------------------##### Collinearity Analysis #####------------------------------------------------------------------
library(caret)
library(randomForest)
library(terra)
library(corrplot)
library(Hmisc)

# load data into data frame
training_fielddata  <- read.csv('training_fielddata_1.csv')
training_fielddata$ID <- seq.int(nrow(training_fielddata))

# create spat Vector data from training data
training_pts_vect <- vect(training_fielddata, crs="EPSG:21037", geom=c("X", "Y"))
testing_pts_vect <- vect(testing_fielddata, crs="EPSG:21037", geom=c("X", "Y"))

### collinearity analysis here exemplified for CV - same done for mean, standard deviation and rao's q
# load and ectract data
cv_rast <- rast("cv_allLayers.tif")
cv_rast_df <- data.frame(cv_rast)
cv_points_extract <- extract(cv_rast, training_pts_vect, raw = T, ID = T)
cv_points_df <- data.frame(cv_points_extract)
cv_points_df <- na.omit(cv_points_df)
cv_points_matrix <- data.matrix(cv_points_df)
# merge with field data
fielddiv_cv <- merge(cv_points_df, training_fielddata, by="ID")
# loop through all layers of the spectral heterogeneity raster and compare the feature importance for each pair of bands
for (i in 1:28) {
  for (j in 1:28) {
    if (i == j) {
      print("same")
    } else {
      # extract bands from data frame
      df_cv_pair <- data.frame(cv_points_df[,i], cv_points_df[,j])
      colnames(df_cv_pair) <- c(i, j)
      print(c(i, j))
      # determine the correlation coefficient between the two bands
      cor_rho <- cor(df_cv_pair[1], df_cv_pair[2], method = "spearman")
      print(cor_rho[1, 1])
      # threshold test for collinearity - values tested: 0.8 
      if (cor_rho[1, 1] >= 0.8) {
        fielddiv_cv_cor <- data.frame(fielddiv_cv$species_richness, df_cv_pair[,1], df_cv_pair[,2])
        # run default random forest model modelling species richness
        rf_cv_pair <- randomForest(fielddiv_cv.species_richness ~ ., data = fielddiv_cv_cor, importance = TRUE)
        feat_imp <- data.frame(importance(rf_cv_pair))
        feat_imp <- feat_imp$X.IncMSE
        # repeat 100 times to overcome random nature of RF
        for (l in 1:99){
          rf_cv_pair <- randomForest(fielddiv_cv.species_richness ~ ., data = fielddiv_cv_cor, importance = TRUE)
          feat_imp_loop <- data.frame(importance(rf_cv_pair))
          feat_imp <- feat_imp + feat_imp_loop$X.IncMSE
        }
        # determine which band has higher average feature importance
        print(feat_imp)
        if (feat_imp[1] > feat_imp[2]){
          cat("Better band of", i, " and ", j, " is: ", i)
        } else {
          cat("Better band of", i, " and ", j, " is: ", j)
        }
      }
    }
  }
}

# all bands were then noted and manually removed based on the results of this analysis

----------------------------------##### Random Forest Models #####---------------------------------------------------------------
library(caret)
library(randomForest)
library(terra)
library(corrplot)
library(Hmisc)
library(stats)
library(doParallel)
library(ggplot2)    # Graphics library
library(sf)         # Spatial data types and handling
library(spdep)      # Diagnosing spatial dependence
library(spatialreg) # Spatial lag and spatial error model
library(gstat)

# enable parallel processin
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

# load data into data frame
training_fielddata  <- read.csv('training_fielddata_1.csv')
training_fielddata$ID <- seq.int(nrow(training_fielddata))

testing_fielddata  <- read.csv('testing_fielddata_1.csv')
testing_fielddata$ID <- seq.int(nrow(testing_fielddata))

# create spat Vector data from training data
training_pts_vect <- vect(training_fielddata, crs="EPSG:21037", geom=c("X", "Y"))
testing_pts_vect <- vect(testing_fielddata, crs="EPSG:21037", geom=c("X", "Y"))

# Read in the data and remove N/A records from training data
##### here again examplified by the mean raster but same approach taken for all other spectral heterogeneity metrics
# load and extract data
mean_rast <- rast("mean_allLayers_treemask1_stretched.tif")
names(mean_rast) <- c("B2", "B3", "B4","B5", "B6",
                    "B7", "B8", "B8A", "B11", "B12",
                    "VV", "VH", "NDVI", "NDWI", "NDPI", "EVI", "EVI2", "NBR", "NBR2",
                    "NIRv", "PI", "Clre", "CCI", "IRECI", "NDI45", "SR", "SAVI")
# mean_rast_df <- data.frame(mean_rast)
mean_points_extract <- extract(mean_rast, training_pts_vect, raw = T, ID = T)
mean_points_df <- data.frame(mean_points_extract)
# clean data of missing values through na.omit
mean_points_df <- na.omit(mean_points_df)
mean_points_matrix <- data.matrix(mean_points_df)
# test data
mean_testing <- extract(mean_rast, testing_pts_vect, raw=T, ID=T)
mean_testing_df <- data.frame(mean_testing)
mean_testing_df <- na.omit(mean_testing_df)
mean_fielddiv_test <- merge(mean_testing_df, testing_fielddata, by="ID")
fielddiv_mean <- merge(mean_points_df, training_fielddata, by="ID")

#to show how the combination models were built
fielddiv_mean <- merge(mean_points_df, training_fielddata, by="ID")
mean_testing <- extract(mean_rast, testing_pts_vect, raw=T, ID=T)
mean_testing_df <- data.frame(mean_testing)
mean_testing_df <- na.omit(mean_testing_df)

mean_fielddiv_test <- merge(mean_testing_df, testing_fielddata, by="ID")
mean_cv_testing_df <- merge(mean_testing_df, cv_testing_df, by="ID")
mean_cv_fielddiv_test <- merge(mean_cv_testing_df, testing_fielddata, by="ID")


##### RANDOM FOREST MODELS #####
# build random forest model: rf <- randomForest(target ~ features, data = , mtry = , ntrees= , importance = TRUE)
# set default settings for recursive feature elimination
# subsets of features to be tested
number_features <- c(1:27)
set.seed(42)
# set controls for RFE - 100 repeats for k-fold cross-validation (method=repeatedmean) (number of folds: number = 5)
ctrl <- rfeControl(functions = rfFuncs, method = 'repeatedcv',
                   repeats = 100, number = 5, verbose = F)

##### mean ##### 
# with species richness and all band - same process done for Shannon index and Inverse Simpson index (Inv_Simpson)
# remove columns that are not relevant 
mean_data_rf <- subset(fielddiv_mean, select = -c(ID, Plot, X, Y, Shannon, Inv_Simpson, nr_trees, residuals, fitted))

## remove features from collinearity analysis
#correlation > 0.8
mean_data_selection <- subset(mean_data_rf, select = c(B3, B12, VH, 
                                                       NDWI, NBR2, Clre, NDI45, SAVI, Shannon))
mean_test_data_selection <- subset(mean_testing_df, select = c(B3, B12, VH, 
                                                               NDWI, NBR2, Clre, NDI45, SAVI))
# only features without response variable
mean_features_sel <- subset(mean_data_selection, select = -c(species_richness))
# Run Recursive feature elimination - rfe function of caret package (x-predictors, y-response variable, sizes-list of number of features tested)
# number_features <- c(1:17)
mean_selected_rfe_profile <- rfe(x = mean_features_sel, 
                               y = mean_data_rf$Shannon,
                               sizes = number_features,
                               rfeControl = ctrl)
# show profile for output
mean_selected_rfe_profile
# plot RFE (RMSE over the number of features)
plot(mean_selected_rfe_profile, main="a) mean")
# determine the remaining predictors/features and variable importance
predictors(mean_selected_rfe_profile)
varImp(mean_selected_rfe_profile)
# determine model performance on test data after removal of bands through RFE
postResample(predict(mean_selected_rfe_profile, mean_testing_df[,2:28]), mean_fielddiv_test$Shannon)

# plot feature importance
varimp_mean_sel <- data.frame(feature = row.names(varImp(mean_selected_rfe_profile))[1:8],
                              importance = varImp(mean_selected_rfe_profile)[1:8, 1])
ggplot(data = varimp_mean_sel, 
       aes(x = reorder(feature, -importance), y = importance, fill = feature)) +
  geom_bar(stat="identity", fill = "darkgreen") + labs(x = "Features (band mean)", y = "Variable Importance (%IncMSE)") + 
  geom_text(aes(label = round(importance, 2)), vjust=1.6, color="white", size=4) + 
  theme_bw() + theme(legend.position = "none" + scale_fill_brewer(palette = "PuBuGn"))


# run the hyperparameter tuning
# create grid for extensive grid search on the mtry 
grid_mean_sel <- expand.grid(
  mtry=seq(1, 8))
ctrl_tune <- trainControl(method = "repeatedcv", number=5, repeats = 100, search = "grid")
set.seed(42)

# loop through number of trees (50 to 1000 in 50 tree increments)
for (t in 1:20) {
  trees <- t*50
  print(trees)
  rf_mean_tuned_sel <- train(species_richness ~ ., data = mean_data_selection, method = "rf",
                             tuneGrid = grid_mean_sel, trControl = ctrl_tune, ntree = trees)
  print(rf_mean_tuned_sel)
}
# this determines the optimal ntree and mtry by lowest RMSE

# test with test data - tuned model is run to predict for test data points
# here repeated 1000 times for sensitivity analysis
set.seed(NULL)
rmse_list <- c()
adj_r_squared_list <- c()
intercept_list <- c()
slope_list <- c()
bias_list <- c()
meandiv_list <- c()
pval_list <- c()
test_r2 <- c()

for (i in 1:1000){
  rf_mean_tuned_sel <- randomForest(species_richness ~ ., data = mean_data_selection, ntree = 400, mtry =1, importance = TRUE)
  # to determine the explained variance
  rf_mean_tuned_sel_2 <- randomForest(species_richness ~ ., data = mean_data_selection, ntree = 400, mtry = 1, importance = TRUE, 
                                      xtest = mean_test_data_selection, ytest = mean_fielddiv_test$species_richness)
  predicted_mean_sel <- predict(rf_mean_tuned_sel, mean_testing_df[,2:28], predict.sel=T)
  meandiv_list[i] <- mean(predicted_mean_sel)
  # Calculate the root mean squared error (RMSE)
  rmse_mean_tuned <- sqrt(mean((predicted_mean_sel - mean_fielddiv_test$species_richness)^2))
  rmse_list[i] <- rmse_mean_tuned
  bias_mean_tuned <- mean((mean_fielddiv_test$species_richness - predicted_mean_sel))
  bias_list[i] <- bias_mean_tuned
  # determine adjusted R², slope and intercept for linear model (lm) comparing predicted and observed species diversity
  lm_mean_test_sel <- lm(predicted_mean_sel ~ mean_fielddiv_test$species_richness)
  adjusted_r <- summary(lm_mean_test_sel)$adj.r.squared 
    adj_r_squared_list[i] <- adjusted_r
  intercept_list[i] <- coef(lm_mean_test_sel)[1]
  slope_list[i] <- coef(lm_mean_test_sel)[2]
  pval_list[i] <- summary(lm_mean_test_sel)$coefficients[2, 4]
  test_r2[i] <- rf_mean_tuned_sel_2$test$rsq[length(rf_mean_tuned_sel_2$test$rsq)]
}
# determine mean and standard deviation on all performance metrics
print("RMSE - mean, stdv")
mean(rmse_list)
sd(rmse_list)
print("Adjusted R² - mean, stdv")
mean(adj_r_squared_list)
sd(adj_r_squared_list)
print("Intercept")
mean(intercept_list)
sd(intercept_list)
print("Slope")
mean(slope_list)
sd(slope_list)
print("Bias")
mean(bias_list)
sd(bias_list)
print("Mean")
mean(meandiv_list)
sd(meandiv_list)
print("explained variance")
mean(test_r2)
sd(test_r2)

# run RF model for predicting the full raster
set.seed(3)
rf_mean_tuned_sel <- randomForest(species_richness ~ ., data = mean_data_selection, ntree = 250, mtry = 1, importance = TRUE)
set.seed(3)
rf_mean_tuned_sel_2 <- randomForest(species_richness ~ ., data = mean_data_selection, ntree = 250, mtry = 1, importance = TRUE, 
                                    xtest = mean_test_data_selection, ytest = mean_fielddiv_test$species_richness)
# determine the performance metrics for this specific model
predicted_mean_sel <- predict(rf_mean_tuned_sel, mean_testing_df[,2:28], predict.sel=T)
predicted_mean_training <- predict(rf_mean_tuned_sel, mean_data_rf[,1:27])

# Calculate the root mean squared error (RMSE)
rmse_mean_tuned <- sqrt(mean((predicted_mean_sel - mean_fielddiv_test$Shannon)^2))

lm_mean_test_sel <- lm(predicted_mean_sel ~ mean_fielddiv_test$Shannon)
summary(lm_mean_test_sel)
# plot linear model comparing predicted and observed species diversity
# Compute shared axis limits
lims <- range(c(mean_fielddiv_test$Shannon, predicted_mean_sel), na.rm = TRUE)
# Extract statistics
intercept <- coef(lm_mean_test_sel)[1]
slope <- coef(lm_mean_test_sel)[2]
r2 <- summary(lm_mean_test_sel)$adj.r.squared
pval <- summary(lm_mean_test_sel)$coefficients[2, 4]  

# Plot with equal axis ranges
plot(mean_fielddiv_test$Shannon, predicted_mean_sel,
     main = "b) Shannon Index",
     xlab = "Observed Shannon",
     ylab = "Predicted Shannon",
     col = "blue", pch = 19,
     xlim = lims, ylim = lims)

abline(lm_mean_test_sel, col = "red", lwd = 2)
abline(0, 1, col="black", lwd=2)
text( x = lims[1] ,
      y = lims[2] - 0.15 * diff(lims),
      labels = paste0("y=(", round(mean(slope_list), 2), "±", round(sd(slope_list), 2), ")x + (", round(mean(intercept_list), 2), "±", round(sd(intercept_list), 2),  ")\n",
                      "Adjusted R² = ", round(mean(adj_r_squared_list), 2), "±", round(sd(adj_r_squared_list), 2), "\n",
                      "RMSE = ", round(mean(rmse_list), 2), " ± ", round(sd(rmse_list), 2), "\n",
                      #                      "RMSE = 6.00", " ± ", round(sd(rmse_list), 2), "\n",
                      "Bias = ", round(mean(bias_list), 2), "±", round(sd(bias_list), 2) ),
      pos = 4, cex = 1)

# predict for full raster
# aggregate raster to 30x30 m cells fitting the field plots 
original_rast <- rast('stacked_sentinel_complete.tif')
names(original_rast) <- c("B2", "B3", "B4","B5", "B6",
                          "B7", "B8", "B8A", "B11", "B12",
                          "VV", "VH", "NDVI", "NDWI", "NDPI", "EVI", "EVI2", "NBR", "NBR2",
                          "NIRv", "PI", "Clre", "CCI", "IRECI", "NDI45", "SR", "SAVI")
aggregated_rast_mean <- aggregate(original_rast, fact = 3)
plot(aggregated_rast_mean[[3]])

# predict species richness with tuned and trained RF forest on aggregated raster
predicted_richness_rast <- predict(object = aggregated_rast_mean, model = rf_mean_tuned_sel, cores=4, cpkgs="randomForest")
plot(predicted_richness_rast)
writeRaster(predicted_richness_rast, "modelled_Shannon_seed3_mean1.tif")
# the same is repeated for all other heterogeneity metrics

### example for combination model
### MEAN + CV ###
number_features <- c(1:54)
mean_cv_df <- merge(mean_points_df, cv_points_df,  by="ID")
fielddiv_mean_cv <- merge(mean_cv_df, training_fielddata, by="ID")
# fielddiv_mean_cv
# here x = mean, y=cv
mean_cv_data_rf <- subset(fielddiv_mean_cv, select = -c(ID, Plot, X, Y, Shannon, 
                                                        Inv_Simpson, nr_trees, residuals, fitted))
mean_cv_features_all <- subset(mean_cv_data_rf, select = -c(species_richness))

rf_mean_cv_default_all <- randomForest(species_richness ~ ., data = mean_cv_data_rf, importance = TRUE)
print(rf_mean_cv_default_all)
# remove bands with collinearity and lower feature importance - mean features
mean_cv_data_selection <- subset(mean_cv_data_rf, select = -c(B2.x, B4.x, B5.x, B6.x, 
                                                              B7.x, B8.x, B8A.x, B11.x, 
                                                              VV.x, NDVI.x, NDPI.x, EVI.x, EVI2.x, NBR.x, PI.x, CCI.x,
                                                              NIRv.x, IRECI.x, SR.x))
# cv features
mean_cv_data_selection <- subset(mean_cv_data_selection, select = -c(B7.y, B8.y, B8A.y, 
                                                                     NDPI.y, EVI.y, EVI2.y, NBR.y, NIRv.y, SR.y, SAVI.y))
mean_cv_features_sel <- subset(mean_cv_data_selection, select = -c(species_richness))

# Recursive feature elimination
mean_cv_selected_rfe_profile <- rfe(x = mean_cv_features_sel, 
                                    y = mean_cv_data_rf$species_richness,
                                    sizes = number_features,
                                    rfeControl = ctrl)
print("RFE for Mean and CV - selected bands")
mean_cv_selected_rfe_profile
plot(mean_cv_selected_rfe_profile, main="e) Mean & CV")
predictors(mean_cv_selected_rfe_profile)
varImp(mean_cv_selected_rfe_profile)
postResample(predict(mean_cv_selected_rfe_profile, mean_cv_testing_df[,2:55]), mean_cv_fielddiv_test$species_richness)

# remove features due to RFE
mean_cv_data_selection_rfe <- subset(mean_cv_data_selection, select = c( VH.x, species_richness))
# hyperparameter tuning
grid_mean_cv_sel <- expand.grid(
  mtry=seq(1, 25))
ctrl_tune <- trainControl(method = "repeatedcv", number=5, repeats = 100, search = "grid")
set.seed(42)
for (t in 1:20) {
  trees <- t*50
  print(trees)
  rf_mean_cv_tuned_sel <- train(species_richness ~ ., data = mean_cv_data_selection_rfe, method = "rf",
                                tuneGrid = grid_mean_cv_sel, trControl = ctrl_tune, ntree = trees)
  print(rf_mean_cv_tuned_sel)
}
# test with test data
rf_mean_cv_tuned_sel <- randomForest(species_richness ~ ., data = mean_cv_data_selection_rfe, ntree = 550, mtry = 1, importance = TRUE)
importance(rf_mean_cv_tuned_sel)
varImpPlot(rf_mean_cv_tuned_sel)
predicted_mean_cv_sel <- predict(rf_mean_cv_tuned_sel, mean_cv_testing_df[,2:55], predict.sel=T)
mean_cv_test_cor_sel <- rcorr(predicted_mean_cv_sel, mean_cv_fielddiv_test$species_richness)
mean_cv_test_cor_sel$P
lm_mean_cv_test_sel <- lm(predicted_mean_cv_sel ~ mean_cv_fielddiv_test$species_richness)
summary(lm_mean_cv_test_sel)

### same approach for spectral species diversity data - this was repeated for each number of clusters (20, 50, 100)
#### Spectral Species Diversity ####
spectral_richness_rast <- rast("richness_mean_20.tiff")
spectral_shannon_rast <- rast("shannon_mean_20.tiff")
spectral_simpson_rast <- rast("simpson_mean_20.tiff")
spectral_species_div <- c(spectral_richness_rast, spectral_shannon_rast, spectral_simpson_rast)

ssdiv_points_extract <- extract(spectral_species_div, training_pts_vect, raw = T, ID = T)
ssdiv_points_extract
ssdiv_points_df <- data.frame(ssdiv_points_extract)
# clean data of missing values through na.omit
ssdiv_points_df <- na.omit(ssdiv_points_df)
ssdiv_points_matrix <- data.matrix(ssdiv_points_df)

fielddiv_ssdiv <- merge(ssdiv_points_df, training_fielddata, by="ID")
mean_ssdiv_testing_df <- merge(mean_testing_df, ssdiv_points_df, by="ID")

ssdiv_testing <- extract(spectral_species_div50, testing_pts_vect, raw=T, ID=T)
ssdiv_testing_df <- data.frame(ssdiv_testing)
ssdiv_testing_df <- na.omit(ssdiv_testing_df)
ssdiv_fielddiv_test <- merge(ssdiv_testing_df, testing_fielddata, by="ID")

ssdiv_features_rf <- subset(fielddiv_ssdiv, select = -c(ID, Plot, X, Y, species_richness, Inv_Simpson, nr_trees, residuals, fitted))
ssdiv_features <- subset(ssdiv_features_rf, select = -c(Shannon)) 
# test of all three bands to determine best predictor
# number_features <- c(1:3)
ssdiv_all_rfe_profile <- rfe(x = ssdiv_features, 
                             y = ssdiv_features_rf$Shannon,
                             sizes = number_features,
                             rfeControl = ctrl)
print("RFE for SSDIV - selected bands")
ssdiv_all_rfe_profile
predictors(ssdiv_all_rfe_profile)
varImp(ssdiv_all_rfe_profile)
# determine performance with test data
postResample(predict(ssdiv_all_rfe_profile, ssdiv_testing_df[,2:4]), ssdiv_fielddiv_test$Shannon)

# run combination with mean model
mean_ssdiv_df <- merge(mean_points_df, ssdiv_points_df,  by="ID")
fielddiv_mean_ssdiv <- merge(mean_ssdiv_df, training_fielddata, by="ID")
mean_ssdiv_fielddiv_test <- merge(mean_ssdiv_testing_df, testing_fielddata, by="ID")
# here x = mean, y=ssdiv
mean_ssdiv_data_rf <- subset(fielddiv_mean_ssdiv, select = -c(ID, Plot, X, Y, species_richness, 
                                                              Inv_Simpson, nr_trees, residuals, fitted))
mean_ssdiv_features_all <- subset(mean_ssdiv_data_rf, select = -c(Shannon))
mean_ssdiv_data_rf
mean_ssdiv_data_selection <- subset(mean_ssdiv_data_rf, select = -c(B2 , B4 , B5 , B6 , 
                                                                    B7 , B8 , B8A , B11 , 
                                                                    VV , NDVI , NDPI , EVI , EVI2 , NBR , PI , CCI ,
                                                                    NIRv , IRECI , SR ))
mean_ssdiv_data_selection <- subset(mean_ssdiv_data_selection, select = -c(mean.simpson, mean.richness))
mean_ssdiv_features_sel <- subset(mean_ssdiv_data_selection, select = -c(Shannon))

# RFE
mean_ssdiv_selected_rfe_profile <- rfe(x = mean_ssdiv_features_sel, 
                                       y = mean_ssdiv_data_rf$Shannon,
                                       sizes = number_features,
                                       rfeControl = ctrl)
print("RFE for Mean and ssdiv - selected bands")
mean_ssdiv_selected_rfe_profile
plot(mean_ssdiv_selected_rfe_profile, main="h) Mean & Spectral Species Shannon")
predictors(mean_ssdiv_selected_rfe_profile)
varImp(mean_ssdiv_selected_rfe_profile)
postResample(predict(mean_ssdiv_selected_rfe_profile, mean_ssdiv_testing_df[,2:31]), mean_ssdiv_fielddiv_test$Shannon)
# hyperparameter tuning if applicable
grid_mean_ssdiv_sel <- expand.grid(
  mtry=seq(1, 8))
ctrl_tune <- trainControl(method = "repeatedcv", number=5, repeats = 100, search = "grid")
set.seed(42)
for (t in 1:20) {
  trees <- t*50
  print(trees)
  rf_mean_ssdiv_tuned_sel <- train(Shannon ~ ., data = mean_ssdiv_data_selection, method = "rf",
                                   tuneGrid = grid_mean_ssdiv_sel, trControl = ctrl_tune, ntrees = trees)
  print(rf_mean_ssdiv_tuned_sel)
}
# test with test data
rf_mean_ssdiv_tuned_sel <- randomForest(Shannon ~ ., data = mean_ssdiv_data_selection, ntrees = 550, mtry = 1, importance = TRUE)
importance(rf_mean_ssdiv_tuned_sel)
varImpPlot(rf_mean_ssdiv_tuned_sel)
predicted_mean_ssdiv_sel <- predict(rf_mean_ssdiv_tuned_sel, mean_ssdiv_testing_df[,2:31], predict.sel=T)
mean_ssdiv_test_cor_sel <- rcorr(predicted_mean_ssdiv_sel, mean_ssdiv_fielddiv_test$Shannon)
mean_ssdiv_test_cor_sel$P
lm_mean_ssdiv_test_sel <- lm(predicted_mean_ssdiv_sel ~ mean_ssdiv_fielddiv_test$Shannon)
summary(lm_mean_ssdiv_test_sel)
