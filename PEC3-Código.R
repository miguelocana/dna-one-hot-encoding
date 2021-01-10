##################################### PARÁMETROS #######################################
# PÁRAMETROS
# 
#
#
#


#######################################################################################

###################################### CÓDIGO #########################################


onehot <- function (df, col_seq, onehot_cols_name='H', remove_col_seq = FALSE) {
  # Se asegura de que todas las seq tengan el mismo tamaño (de lo contrario salta ERROR)
  f <- length(strsplit(df[col_seq][[1]], "")[[1]])
  for (j in 1:nrow(df)) {
    a = length(strsplit(df[j,col_seq], "")[[1]])
    if (a == f) {
      next
    } else {
      print('ERROR: al menos una secuencia no tiene el mismo tamaño que el resto.')
      break
    }
  }
  
  # Crea una columna por cada dígito del one-hot
  num_cols <- f * 4
  cols_onehot <- list()
  for (i in 1:num_cols) {
    cols_onehot <- append(cols_onehot,paste(onehot_cols_name,i,sep=""))
  }
  
  # Añade las columnas (vacías) al nuevo DF
  df2 <- df
  for (i in cols_onehot) {
    df2[,i] <- NA
  }
  
  # Traducción de seq a one-hot encoding
  for (j in 1:nrow(df2)) {
    seq = strsplit(tolower(df2[j,col_seq]), "")[[1]] # convierte la seq en minuscula
    seq_onehot <- list()
    for (i in seq) {
      if (i == 'a') {
        seq_onehot <- append(seq_onehot, c(0,0,0,1))
      } else if (i == 'g') {
        seq_onehot <- append(seq_onehot, c(0,0,1,0))
      } else if (i == 'c') {
        seq_onehot <- append(seq_onehot, c(0,1,0,0))
      } else if (i == 't') {
        seq_onehot <- append(seq_onehot, c(1,0,0,0))
      } else {
        print("ERROR: al menos un carácter de la secuencia introducida no corresponde a ningún nucleótido (A,G,C,T).")
        break
      }
    }
    for (k in 1:length(seq_onehot)) {
      ind <- grep(paste(onehot_cols_name,k,sep=""), colnames(df2))
      df2[j,ind] <- seq_onehot[k]
    }
  }
  
  # Elimina o no la columna con la secuencia
  if (remove_col_seq) {
    df3 <- df2[,c(x[-grep(col_seq, colnames(df2))])]
  } else {
    df3 <- df2
  }
}
