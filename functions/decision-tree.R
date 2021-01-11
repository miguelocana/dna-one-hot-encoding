##################################### LIBRERÍAS ########################################
library('C50') # Decision Tree

# ARBOL DE DECISION, observaciones:

# Carga datos
promoters <- read.csv('promoters.txt',header = FALSE,sep=',')
df <- onehot(df=promoters[,-2], col_seq = 'V3', onehot_cols_name = 'P.', remove_col_seq = TRUE)

