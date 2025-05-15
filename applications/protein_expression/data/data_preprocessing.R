# Preprocessing Data as in SOM paper: (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0129126)

library(readxl)
library(stringr)
library(dplyr)
library(naniar)
library(ggplot2)
library(tidyr)
library(readr)

#
# read in raw data -------------
#

df <- data.frame(read_xls("Data_Cortex_Nuclear.xls"))

# split up mouse id and trial number

df[c("Mouse_ID", "Trial")] <- str_split_fixed(df$MouseID, '_', 2) 

df$MouseID <- NULL
df <- df[ , c(82,83,1:81)]

###############################################################
#
# look for outlier mouse -------------------------------------------------------------------
#

# total missingness
sum(is.na(df))  # 1396

# protein missingness
miss_var_summary(df)
temp <- colSums(is.na(df[ ,c(3:79)]))
temp
table(temp)

# mouse missingness
temp <- rowSums(is.na(df[ ,c(3:79)]))
temp
table(temp)

# potential for outlier from missingness - index 988-990 (mouse 3426)

# compare protein expression values between classes

length(unique(df$Mouse_ID))  # 72
table(df$class)/15


t <- df %>%
  filter(class == "t-CS-m")

# plot by protein

ggplot(t, aes(x=Mouse_ID, y=pAKT_N)) +
  geom_boxplot()


t_long <- t %>%
  select(c(1,2,63:72)) %>%
  pivot_longer(cols = 3:12, names_to = "protein", values_to = "value")


ggplot(t_long, aes(x=Mouse_ID, y=value)) +
  geom_boxplot() + 
  facet_grid(rows = vars(protein), scales = "free")


# Conclusion: 
#  - could not find the outlier mouse mentioned in the paper.  Maybe this mouse was already removed
#    from the dataset?
#  - SOM uses all 72 mice (and 15 measurements per mouse)
#  - Will NOT remove any mice from the set of 72
###########################################################################################



#
# fill in missing values in data -----------------------------
#

# replace missing values with average value from class

df_fill <- data.frame(df)

classes <- unique(df_fill$class)
offset = 0

for(c in classes)
{
  # get row indices for class c
  row_idx <- which(df_fill$class == c)
  
  # loop through columns (proteins)
  for(i in 3:79)
  {
    na_idx = offset + which(is.na(df_fill[row_idx, i]))
    df_fill[na_idx, i] <- mean(df_fill[row_idx, i], na.rm = TRUE)
  }
  
  offset = row_idx[length(row_idx)]
}


#
# min-max normalization ---------------------------
#

df_fill_norm <- data.frame(df_fill)

# loop through each column (protein)

for(i in 3:79)
{
  col_min <- min(df_fill_norm[ , i])
  col_max <- max(df_fill_norm[ , i])
    
  df_fill_norm[ , i] <- (df_fill_norm[ , i] - col_min)/(col_max - col_min)   
  
}


#
# write preprocessed data ---------------------------------
#

write_csv(df_fill_norm, "Preprocessed_Data.csv")
