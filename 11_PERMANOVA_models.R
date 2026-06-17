## -----------------------------------------------------------------------------#
# Project: "Taxonomic and functional composition of mycorrhizal communities respond 
# differently to host identity and environment"
#
# New PERMANOVA models using environmental PCA axis values 
# 
# Original author: L. McKinley Nevins 
# 
# April 29, 2026
#
# Software versions:  R v 4.5.2
#                     tidyverse v 2.0.0
#                     vegan v 2.7.3
#                     dplyr v 1.1.4
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")

#################################################################################
#                               Main workflow                                   #
#  Create new PERMANOVA models to assess the effects of host identity and site  #
#  environmental factors on AM and ECM taxonomic and functional dissimilarity,  #
#  using the new environmental PCA axis values for each host tree fungal        #
#  community.                                                                   #
#################################################################################


################ --
# (1) DATA PREP
################ --

wd <- "~/Dropbox/WSU/Mycorrhizae_Project/"
setwd(wd)


# Load in Aitchison distance matrices for both AM and ECM communities 

# For taxa
load("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_AM_2025.Rdata")

AM_tax_mat <- as.matrix(aitchison_AM)


load("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_EM_no_THPL.Rdata")

ECM_tax_mat <- as.matrix(aitchison_EM)


# For functions
load("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_AM_traits.Rdata")

AM_func_mat <- as.matrix(aitchison_AM_traits)


load("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/aitchison_dist_EM_traits_no_THPL_2.0.Rdata")

ECM_func_mat <- as.matrix(aitchison_EM_traits)


# Load in environmental PC1, PC2, and PC3 values for each Field_ID site 
all_enviro_PCs <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/PCA_Outputs/all_enviro_PCs.csv")


# Load in full environmental datasets for both to get sample codes 
AM_env <- read.csv(file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_enviro_all_2025.csv")

ECM_env <- read.csv(file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_enviro_no_THPL.csv")


# Subset both 
AM_env <- dplyr::select(AM_env, Sample_Code = X, Field_ID, Host_ID)

ECM_env <- dplyr::select(ECM_env, Sample_Code = X, Field_ID, Host_ID)


# Merge to PC dataframes 
AMtrees_enviro_PCs <- merge(all_enviro_PCs, AM_env, by = "Field_ID")

EMtrees_enviro_PCs <- merge(all_enviro_PCs, ECM_env, by = "Field_ID")



# Set Sample_Code to the rownames for both 
AMtrees_enviro_PCs <- column_to_rownames(AMtrees_enviro_PCs, "Sample_Code")

EMtrees_enviro_PCs <- column_to_rownames(EMtrees_enviro_PCs, "Sample_Code")

#################################################################################

######################### --
# (2) PERMANOVA MODELING
######################### --


### AM ###

## TAXONOMIC ## 

permanova.AM_tax_enviro <- vegan::adonis2(aitchison_AM ~ PC1 + PC2 + PC3, data = AMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.AM_tax_enviro


#           Df SumOfSqs      R2      F Pr(>F)    
# Model      3    62025 0.09492 4.4045  0.001 ***
# Residual 126   591448 0.90508                  
# Total    129   653472 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Explained 9.5% 


permanova.AM_tax_host <- vegan::adonis2(aitchison_AM ~ Host_ID, data = AMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.AM_tax_host


#           Df SumOfSqs     R2      F Pr(>F)    
# Model      2    22088 0.0338 2.2215  0.001 ***
# Residual 127   631384 0.9662                  
# Total    129   653472 1.0000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Explained 3.4%


permanova.AM_tax_both <- vegan::adonis2(aitchison_AM ~ PC1 + PC2 + PC3 + Host_ID, data = AMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.AM_tax_both

#           Df SumOfSqs      R2      F Pr(>F)    
# Model      5    84731 0.12966 3.6947  0.001 ***
# Residual 124   568742 0.87034                  
# Total    129   653472 1.00000               
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# explained 13.0%


permanova.AM_tax_both2 <- vegan::adonis2(aitchison_AM ~ PC1 + PC2 + PC3 + Host_ID + PC1*Host_ID + PC2*Host_ID + PC3*Host_ID, data = AMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.AM_tax_both2

#           Df SumOfSqs      R2      F Pr(>F)    
# Model     11   140074 0.21435 2.9268  0.001 ***
# Residual 118   513399 0.78565                  
# Total    129   653472 1.00000                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# explained 21.4%



## FUNCTIONAL ## 

permanova.AM_func_enviro <- vegan::adonis2(aitchison_AM_traits ~ PC1 + PC2 + PC3, data = AMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.AM_func_enviro


#           Df SumOfSqs      R2      F Pr(>F)
# Model      3    93.27 0.03181 1.3797  0.219
# Residual 126  2839.12 0.96819              
# Total    129  2932.39 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Explained 3.2% 


permanova.AM_func_host <- vegan::adonis2(aitchison_AM_traits ~ Host_ID, data = AMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.AM_func_host


#           Df SumOfSqs      R2    F Pr(>F)   
# Model      2    187.6 0.06397 4.34  0.003 **
# Residual 127   2744.8 0.93603               
# Total    129   2932.4 1.00000                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Explained 6.4% 


permanova.AM_func_both <- vegan::adonis2(aitchison_AM_traits ~ PC1 + PC2 + PC3 + Host_ID, data = AMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.AM_func_both

#             Df SumOfSqs      R2      F Pr(>F)   
# Model      5   323.13 0.11019 3.0712  0.002 **
# Residual 124  2609.26 0.88981                 
# Total    129  2932.39 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# explained 11.0%


permanova.AM_func_both2 <- vegan::adonis2(aitchison_AM_traits ~ PC1 + PC2 + PC3 + Host_ID + PC1*Host_ID + PC2*Host_ID + PC3*Host_ID, data = AMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.AM_func_both2

#           Df SumOfSqs      R2     F Pr(>F)   
# Model     11   530.08 0.18077 2.367  0.003 **
# Residual 118  2402.31 0.81923                
# Total    129  2932.39 1.00000                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# explained 18.1%


############################## -- 

### ECM ###

## TAXONOMIC ## 

permanova.ECM_tax_enviro <- vegan::adonis2(aitchison_EM ~ PC1 + PC2 + PC3, data = EMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.ECM_tax_enviro


#           Df SumOfSqs      R2     F Pr(>F)    
# Model      3    26709 0.02061 2.574  0.001 ***
# Residual 367  1269395 0.97939                 
# Total    370  1296104 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Explained 2.1% 


permanova.ECM_tax_host <- vegan::adonis2(aitchison_EM ~ Host_ID, data = EMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.ECM_tax_host


#          Df SumOfSqs      R2      F Pr(>F)    
# Model      6    26259 0.02026 1.2545  0.001 ***
# Residual 364  1269844 0.97974                  
# Total    370  1296104 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Explained 2.0% 


permanova.ECM_tax_both <- vegan::adonis2(aitchison_EM ~ PC1 + PC2 + PC3 + Host_ID, data = EMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.ECM_tax_both

#           Df SumOfSqs      R2      F Pr(>F)    
# Model      9    54323 0.04191 1.7547  0.001 ***
# Residual 361  1241781 0.95809                  
# Total    370  1296104 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# explained 4.2%


permanova.ECM_tax_both2 <- vegan::adonis2(aitchison_EM ~ PC1 + PC2 + PC3 + Host_ID + PC1*Host_ID + PC2*Host_ID + PC3*Host_ID, data = EMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.ECM_tax_both2

#           Df SumOfSqs      R2      F Pr(>F)    
# Model     27   140568 0.10845 1.5454  0.001 ***
# Residual 343  1155535 0.89155                  
# Total    370  1296104 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# explained 10.8%



## FUNCTIONAL ## 

permanova.ECM_func_enviro <- vegan::adonis2(aitchison_EM_traits ~ PC1 + PC2 + PC3, data = EMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.ECM_func_enviro


#          Df SumOfSqs      R2      F Pr(>F)    
# Model      3     2095 0.04027 5.1333  0.001 ***
# Residual 367    49921 0.95973                  
# Total    370    52016 1.00000                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Explained 4.0% 


permanova.ECM_func_host <- vegan::adonis2(aitchison_EM_traits ~ Host_ID, data = EMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.ECM_func_host


#           Df SumOfSqs     R2      F Pr(>F)  
# Model      6     1051 0.0202 1.2506  0.096 .
# Residual 364    50965 0.9798                
# Total    370    52016 1.0000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Explained 2.0% 


permanova.ECM_func_both <- vegan::adonis2(aitchison_EM_traits ~ PC1 + PC2 + PC3 + Host_ID, data = EMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.ECM_func_both

#           Df SumOfSqs      R2      F Pr(>F)    
# Model      9     3280 0.06306 2.6998  0.001 ***
# Residual 361    48736 0.93694                  
# Total    370    52016 1.00000                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# explained 6.3%


permanova.ECM_func_both2 <- vegan::adonis2(aitchison_EM_traits ~ PC1 + PC2 + PC3 + Host_ID + PC1*Host_ID + PC2*Host_ID + PC3*Host_ID, data = EMtrees_enviro_PCs, method = "euclidean", permutations = 999)
permanova.ECM_func_both2

#           Df SumOfSqs      R2      F Pr(>F)    
# Model     27     8179 0.15723 2.3701  0.001 ***
# Residual 343    43837 0.84277                  
# Total    370    52016 1.00000                          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# explained 15.7%


# Interpretation: The power of these interaction terms is essentially showing that the effect of host identity on the 
# community composition depends on the environmental context


## -- END -- ## 
