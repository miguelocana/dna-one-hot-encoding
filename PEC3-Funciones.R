##################################### LIBRER?AS ########################################
library('class') # kNN
library('kernlab') # Support Vector Machine (SVM)
library('e1071') # Naive Bayes
library('randomForest') # Random Forest
library('C50') # Decision Tree
library('neuralnet') # Artifical Neural Network
library('gmodels') 
library('caret') # Confusion Matrix

##################################### PAR?METROS #######################################
promoters <- read.csv('promoters.txt',header = FALSE,sep=',')

#!! - La clase o variable dependiente debe situarse en la primera columna.
#!! - La secuencia de ADN debe situarse en la ?ltima columna.

promoters <- promoters[,-2]

##################################### DATAFRAMES #######################################
df_onehot <- onehot(df = promoters, col_seq = 'V3')
df_categorico <- sepseq(df = promoters, col_seq = 'V3')

###################################### C?DIGO ##########################################
onehot <- function (df, col_seq, onehot_cols_name='H', remove_col_seq = TRUE) {
  # Se asegura de que todas las seq tengan el mismo tama?o (de lo contrario salta ERROR)
  f <- length(strsplit(df[col_seq][[1]], "")[[1]])
  for (j in 1:nrow(df)) {
    a = length(strsplit(df[j,col_seq], "")[[1]])
    if (a == f) {
      next
    } else {
      print('ERROR: al menos una secuencia no tiene el mismo tama?o que el resto.')
      break
    }
  }
  
  # Crea una columna por cada d?gito del one-hot
  num_cols <- f * 4
  cols_onehot <- list()
  for (i in 1:num_cols) {
    cols_onehot <- append(cols_onehot,paste(onehot_cols_name,i,sep=""))
  }
  
  # A?ade las columnas (vac?as) al nuevo DF
  df2 <- df
  for (i in cols_onehot) {
    df2[,i] <- NA
  }
  
  # Traducci?n de seq a one-hot encoding
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
        print("ERROR: al menos un car?cter de la secuencia introducida no corresponde a ning?n nucle?tido (A,G,C,T).")
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
sepseq <- function (df, col_seq, sepseq_cols_name='N', remove_col_seq = TRUE) {
  # Se asegura de que todas las seq tengan el mismo tama?o (de lo contrario salta ERROR)
  f <- length(strsplit(df[col_seq][[1]], "")[[1]])
  for (j in 1:nrow(df)) {
    a = length(strsplit(df[j,col_seq], "")[[1]])
    if (a == f) {
      next
    } else {
      print('ERROR: al menos una secuencia no tiene el mismo tama?o que el resto.')
      break
    }
  }
  
  # Crea una columna por cada nucleotido
  num_cols <- f 
  cols_sepseq <- list()
  for (i in 1:num_cols) {
    cols_sepseq <- append(cols_sepseq,paste(sepseq_cols_name,i,sep=""))
  }
  
  # A?ade las columnas (vac?as) al nuevo DF
  df2 <- df
  for (i in cols_sepseq) {
    df2[,i] <- NA
  }
  
  # Mueve cada nucle?tido a su correspondiente columna
  for (j in 1:nrow(df2)) {
    seq = strsplit(tolower(df2[j,col_seq]), "")[[1]] # convierte la seq en minuscula
    seq_sepseq <- list()
    for (i in seq) {
      if ((i == 'a') | (i == 'g') | (i == 'c') | (i == 't')) {
        seq_sepseq <- append(seq_sepseq, i)
      } else {
        print("ERROR: al menos un car?cter de la secuencia introducida no corresponde a ning?n nucle?tido (A,G,C,T).")
        break
      }
    }
    for (k in 1:length(seq_sepseq)) {
      ind <- grep(paste(sepseq_cols_name,k,sep=""), colnames(df2))
      df2[j,ind] <- seq_sepseq[k]
    }
  }
  
  # Elimina o no la columna con la secuencia
  if (remove_col_seq) {
    df3 <- df2[,c(colnames(df2)[-grep(col_seq, colnames(df2))])]
  } else {
    df3 <- df2
  }
  
  # Devuelve el df con los nucle?tidos separados en columnas diferentes
  return (df3)
}
split_train_test <- function(df, size = 0.7, seed = 123) {
  # Establece la semilla
  set.seed(seed)
  
  # Devuelve un "diccionario" con el DF partido en train y test
  t <- floor(size * nrow(df))
  train_ind <- sample(seq_len(nrow(df)), size = t)
  train <- df[train_ind,]
  test <- df[-train_ind,]
  
  # Quita la semilla
  set.seed(Sys.time())
  
  return (list('train' = train, 'test' = test))
}

##################################### ALGORITMOS ######################################
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
naive_bayes <- function (df_categorico) {
  
  # - Se explorar? la opci?n de activar o no 'laplace'.
  # - Tanto las variables dependientes como independientes deben ser de tipo factor.
  # - Se separan los clasificadores de train y test.
  
  # Transforma todas las columnas en factores
  for (i in 1:ncol(df)) {
    # Comprueba si es factor, de lo contrario lo transforma
    if (is.factor(df[,i])) {
      next
    } else {
      df[,i] <- as.factor(df[,i])
    }
  }
  
  # Split train y data
  train <- split_train_test(df)$train[,-1]
  test <- split_train_test(df)$test[,-1]
  train_labels <- as.factor(split_train_test(df)$train[,1])
  test_labels <- as.factor(split_train_test(df)$test[,1])
  
  # Laplace = 0
  # Modelo y estad?sticos
  NBlp0 <- e1071::naiveBayes(train, train_labels, laplace = 0)
  # Predicciones
  pNBlp0 <- predict(NBlp0, test)
  c1 <- confusionMatrix(pNBlp0, test_labels, dnn = c('Predicho','Actual'))
  # Interpretaci?n final:
  # - El modelo tiene una precisi?n del 90%
  # - Devuelve 3 falsos negativos
  
  # Laplace = 1
  # Modelo y estad?sticos
  NBlp1 <- e1071::naiveBayes(train, train_labels, laplace = 1)
  # Predicciones
  pNBlp1 <- predict(NBlp1, test)
  c2 <- confusionMatrix(pNBlp1, test_labels, dnn = c('Predicho','Actual'))
  # Interpretaci?n final:
  # - El modelo tiene una precisi?n del 84%
  # - Devuelve 2 falsos positivos 
  # - Devuelve 1 falso negativo
  
  return (list('NBlp0'=c1, 'NBlp1'=c2))
}
decision_tree <- function(df_onehot) {
  
  # - Hay que separar la variable dependiente de los predictores en train y test
  # - El vector con la variable dependiente debe ser un factor (al menos con C5.0)
  
  # Split train y data
  train <- split_train_test(df_onehot)$train[,-1]
  test <- split_train_test(df_onehot)$test[,-1]
  train_labels <- as.factor(split_train_test(df_onehot)$train[,1])
  test_labels <- as.factor(split_train_test(df_onehot)$test[,1])
  
  # Sin boosting
  # Modelo y estad?sticos
  DTsb <- C5.0(train, train_labels)
  # Predicciones
  pDTsb <- predict(DTsb, test)
  c1 <- confusionMatrix(pDTsb, test_labels, dnn = c('Predicho','Actual'))
  
  # Con boosting
  # Modelo y estad?sticos
  DTnb <- C5.0(train, train_labels, trials = 10)
  # Predicciones
  pDTnb <- predict(DTnb, test)
  c2 <- confusionMatrix(pDTnb, test_labels, dnn = c('Predicho','Actual'))
  
  return (list('DTsb'=c1, 'DTnb'=c2))
}
random_forest <- function(df_onehot) {
  
  # - Se explorar?n la opci?n de n?mero de ?rboles n - 50, 100
  # - Hay que separar la variable dependiente de los predictores en train y test
  # - La variable dependiente debe ser un factor
  
  # Convierte la clase en factor
  df_onehot$V1 <- as.factor(df_onehot$V1)
  
  # Split train y data y labels
  train <- split_train_test(df_onehot)$train[,-1]
  test <- split_train_test(df_onehot)$test[,-1]
  train_labels <- split_train_test(df_onehot)$train[,1]
  test_labels <- split_train_test(df_onehot)$test[,1]
  
  # RF 50
  # Modelo y estad?sticos
  RF50 <- randomForest(train, train_labels, ntree = 50)
  # Predicciones
  pRF50 <- predict(RF50, test)
  confusionMatrix(pRF50, test_labels, dnn = c('Predicho','Actual'))
  # Interpretaci?n final:
  # - El modelo tiene una precisi?n del 93%
  # - Devuelve 1 falso positivo 
  # - Devuelve 1 falso negativo
  
  # RF 100
  # Modelo y estad?sticos
  RF100 <- randomForest(train, train_labels, ntree = 100)
  # Predicciones
  pRF100 <- predict(RF100, test)
  c2 <- confusionMatrix(pRF100, test_labels, dnn = c('Predicho','Actual'))
  # Interpretaci?n final:
  # - El modelo tiene una precisi?n del 87%
  # - Devuelve 2 falsos positivos
  # - Devuelve 2 falsos negativos
  
  return (list('RF50'=c1, 'RF100'=c2))
}
svm <- function(df_onehot) {
  
  # - Se explorar?n las funciones kernel lineal y rbf.
  # - Las variables independientes deben ser num?ricas y estar normalizadas.
  # - Las variables dependientes deben ser de tipo factor.
  # - Se separan los clasificadores de train y test.
  
  # Transforma la clase en tipo factor
  df_onehot[,1] <- as.factor(df_onehot[,1])
  
  # Split train y data
  train <- split_train_test(df_onehot)$train
  test <- split_train_test(df_onehot)$test
  
  # Kernel lineal
  # Modelo y estad?sticos
  SVM_lineal <- ksvm(V1 ~ ., data = train, kernel = 'vanilladot')
  # Predicciones
  pSVM_lineal <- predict(SVM_lineal, test)
  c1 <- confusionMatrix(pSVM_lineal, test$V1, dnn = c('Predicho','Actual'))
  
  # Intepretaci?n:
  # - El modelo tiene una precisi?n del 90%
  # - Devuelve 2 falsos positivos 
  # - Devuelve 1 falsos negativos
  
  # RBF
  # Modelo y estad?sticos
  SVM_rbf <- ksvm(V1 ~ ., data = train, kernel = 'rbfdot')
  # Predicciones
  pSVM_rbf <- predict(SVM_rbf, test)
  c2 <- confusionMatrix(pSVM_rbf, test$V1, dnn = c('Predicho','Actual'))
  # Intepretaci?n:
  # - El modelo tiene una precisi?n del 90%
  # - Devuelve 1 falso positivo
  # - Devuelve 2 falsos negativos 
  
  return (list('SVM_lineal'=c1, 'SVM_rbf'=c2))
}
ann <- function(df_onehot) {
  
  # - Se explorar?n el n?mero de nodos de la capa oculta n - 4,5.
  # - Se necesita una columna por cada tipo de clasificador.
  # - Las variables dependientes deben ser de tipo num?rico (entre 0 y 1).
  
  # Transforma la clase en tipo factor
  df_onehot[,1] <- as.factor(df_onehot[,1])
  
  # Split train y data
  train <- split_train_test(df_onehot)$train
  test <- split_train_test(df_onehot)$test
  
  # ANN 4
  # Modelo y estad?sticos
  ANN4 <- neuralnet(V1 + ~., data = train_ann,hidden = 4,linear.output = FALSE )
  # Predicciones
  pANN4 <- compute(ANN4,test)$net.result
  c1 <- confusionMatrix(pANN4, test$V1, dnn = c('Predicho','Actual'))
  
  # Intepretaci?n:
  # - El modelo tiene una precisi?n del 90%
  # - Devuelve 2 falsos positivos 
  # - Devuelve 1 falsos negativos
  
  # ANN 5
  # Modelo y estad?sticos
  SVM_rbf <- ksvm(V1 ~ ., data = train, kernel = 'rbfdot')
  # Predicciones
  pSVM_rbf <- predict(SVM_rbf, test)
  c2 <- confusionMatrix(pSVM_rbf, test$V1, dnn = c('Predicho','Actual'))
  # Intepretaci?n:
  # - El modelo tiene una precisi?n del 90%
  # - Devuelve 1 falso positivo
  # - Devuelve 2 falsos negativos 
  
  return (list('SVM_lineal'=c1, 'SVM_rbf'=c2))
}

cat_to_num <- function (df) {
  # Está pensada para la variable dependiente,
  # la primera columna debe ser la clase.
  # Transforma los promotores en 1 (+) y 0 (-)
  for (i in 1:nrow(df)) {
    if (df[i,1]=='+') {
      df[i,1] <- 1
    } else {
      df[i,1] <- 0
    }
  }
  df[,1] <- as.numeric(df[,1])
  return (df)
}
num_to_cat <- function (df) {
  # Está pensada para la variable dependiente,
  # la primera columna debe ser la clase.
  # Transforma los promotores en factores + (1) y - (0)
  for (i in 1:nrow(df)) {
    if (df[i,1]==1) {
      df[i,1] <- '+'
    } else {
      df[i,1] <- '-'
    }
  }
  df[,1] <- as.factor(df[,1])
  return (df)
}

df_onehot <- onehot(df = promoters, col_seq = 'V3')
df_ann <- cat_to_num(df_onehot)
# Split train y data
train <- split_train_test(df_ann)$train
test <- split_train_test(df_ann)$test

# ANN 4
# Modelo y estad?sticos
ANN4 <- neuralnet(V1 ~., data = train,hidden = 4,linear.output = FALSE )
# Predicciones
pANN4 <- round(compute(ANN4,test)$net.result,digits=0)
pANN4

test <- num_to_cat(test)
pANN4 <- num_to_cat(pANN4)
pANN4
c1 <- confusionMatrix(as.factor(round(pANN4, digits = 0)), as.factor(test$V1), dnn = c('Predicho','Actual'))

# Intepretaci?n:
# - El modelo tiene una precisi?n del 84%
# - Devuelve 4 falsos positivos 
# - Devuelve 1 falsos negativos

# ANN 5
# Modelo y estad?sticos
SVM_rbf <- ksvm(V1 ~ ., data = train, kernel = 'rbfdot')
# Predicciones
pSVM_rbf <- predict(SVM_rbf, test)
c2 <- confusionMatrix(pSVM_rbf, test$V1, dnn = c('Predicho','Actual'))
# Intepretaci?n:
# - El modelo tiene una precisi?n del 90%
# - Devuelve 1 falso positivo
# - Devuelve 2 falsos negativos 