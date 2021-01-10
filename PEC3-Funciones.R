##################################### LIBRERÍAS ########################################
library('class')
library('kernlab')
library('e1071')
library('randomForest')
library('C50')
library('neuralnet')
library('gmodels')

##################################### PARÁMETROS #######################################
promoters <- read.csv('promoters.txt',header = FALSE,sep=',')

# !! - La clase siempre debe ser la primera columna
# !! - La secuencia siempre debe ser la última columna

###################################### CÓDIGO ##########################################
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
    df3 <- df2[,c(colnames(df2)[-grep(col_seq, colnames(df2))])]
  } else {
    df3 <- df2
  }
  
  # Devuelve el df con one-hot encoding
  return (df3)
}

split_train_test <- function(df, size = 0.7) {
  # Devuelve un "diccionario" con el DF partido en train y test
  t <- floor(size * nrow(df))
  train_ind <- sample(seq_len(nrow(df)), size = t)
  train <- df[train_ind,]
  test <- df[-train_ind,]
  return (list('train' = train, 'test' = test))
}

kneighbors <- function(train, test) {
  train_labels <- train[,1]
  test_labels <- test[,1]
  kneighbors_1 <- knn(train = train[,-1], test = test[,-1], cl = train_labels, k = 1)
  kneighbors_3 <- knn(train = train[,-1], test = test[,-1], cl = train_labels, k = 3)
  kneighbors_5 <- knn(train = train[,-1], test = test[,-1], cl = train_labels, k = 5)
  kneighbors_7 <- knn(train = train[,-1], test = test[,-1], cl = train_labels, k = 7)
  kneighbors_1_p <- CrossTable(x = test_labels, y = kneighbors_1, prop.chisq=FALSE)
  kneighbors_3_p <- CrossTable(x = test_labels, y = kneighbors_3, prop.chisq=FALSE)
  kneighbors_5_p <- CrossTable(x = test_labels, y = kneighbors_5, prop.chisq=FALSE)
  kneighbors_7_p <- CrossTable(x = test_labels, y = kneighbors_7, prop.chisq=FALSE)
  return (list('k1'=kneighbors_1_p,'k3'=kneighbors_3_p,'k5'=kneighbors_5_p,'k7'=kneighbors_7_p))
}

ann <- function(train, test) {
  
}

##################################### ALGORITMOS ######################################
ann <- neuralnet(formula = V1 ~ ., data = train, hidden = 1)
p <- sapply(compute(ann, test)$net.result, round, digits = 0)


svm <- ksvm(target ~predictors, data = mydata, kernal = "rbfdot", c = 1)
p <- predict(m, test, type = 'response')
p$net.result

decision_tree <- c5.0(train, class, trials = 1, costs, NULL)
p <- predict(m, test, type = 'class')

random_forest <- randomForest(train, class, ntree = 500, mtry = sqrt(p))
p <- predict(m, test, type = 'response')

###################################### PRE-PROCESAMIENTO ##############################
df1 <- promoters[,-2] # quitamos la columna "name"
df2 <- onehot(df=df1, col_seq = 'V3', onehot_cols_name = 'P.', remove_col_seq = TRUE)
df2$V1 <- factor(df2$V1, levels = c("+", "-"), labels = c("Promotor", "No promotor"))

###################################### SPLIT ##########################################
train <- split_train_test(prueba)$train
test <- split_train_test(prueba)$test
train_labels <- df_train[,1]
test_labels <- df_test[,1]

##################################### RESULTADOS ####################################
KN <- kneighbors(train, test)$k1
NB <- 

#########################################################################
for (i in 1:ncol(df2)) {
  df2[,i] <- as.character(df2[,i])
}

naive_bayes <- naiveBayes(V1 ~ ., data = train)
pred <- predict(naive_bayes, test)
naive_bayes_cr <- CrossTable(x = naive_bayes_p, y = test_labels, prop.chisq = FALSE)

#########################################################################
