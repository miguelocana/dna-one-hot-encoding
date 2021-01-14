Los promotores son secuencias de ADN que afectan la frecuencia y ubicación del inicio de la transcripción através de la interacción con la ARN polimerasa. El estudio se basa en ficheros obtenidos de:

Dua, D. and Graff, C. (2019). UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. Irvine, CA:University of California, School of Information and Computer Science.

Los atributos del fichero de datos (***promoters.txt***) son:

1. Un símbolo de {+/-}, indicando la clase (“+” = promotor).
2. El nombre de la secuencia promotora. Las instancias que corresponden a no promotores se denominanpor la posición genómica.
3. Las restantes 57 posiciones corresponden a la secuencia

La manera elegida para representar los datos es un paso crucial en los algoritmos de clasificación. En el caso que nos ocupa, análisis basados en secuencias, se usará la transformación denominada *one-hot encoding*.

El ***one-hot encoding*** representa cada nucleótido por un vector de 4 componentes, con 3 de ellas a 0 y una a 1 indicando el nucleótido. Pongamos por ejemplo, el nucleótido T se representa por (1,0,0,0), el nucleótido C por (0,1,0,0), el nucleótido G por (0,0,1,0) y el nucleótido A por (0,0,0,1). Por tanto, para una secuencia de 57 nucleótidos, como en nuestro caso, se obtendrá un vector de 4*57 = 228 componentes, resultado de concatenar los vectores para cada uno de los 57 nucleótidos. 

Una vez realizada la transformación *one-hot encoding*, el objetivo se trata de implementar distintos algoritmos para predecir si la secuencia es un promotor o no, y comparar sus rendimientos.

**Ficheros**:
- ```functions.R```: se encuentran todas las funciones del proyecto.
- ```informe.Rmd```: genera informes dinámicos.

**TODO**: 
- app.R: contiene código para ```shiny```
