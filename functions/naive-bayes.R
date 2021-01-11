##################################### LIBRERÍAS ########################################
library('e1071') # Naive Bayes

# NAIVE BAYES, observaciones:


# Carga datos
promoters <- read.csv('promoters.txt',header = FALSE,sep=',')
df_categorica <- sepseq(df=promoters[,-2], col_seq = 'V3', remove_col_seq = TRUE)


