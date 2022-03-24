
## Microarreglos: expresion diferencial y anotacion

library(Biobase) ## Cargar libreria 
library(limma) ## Cargar libreria 


## Introduccion

expre <- read.table ("expr_normalizada.txt", header = TRUE, row.names = 1) ## Cargar el archivo de texto de expresion
head (expre) ## Ver la primer parte de la tabla

## Para algo no tan funde ojos, el paquete Viridis tiene distintos ajustes para gradientes de color

install.packages ("viridis") ## Instalar el paquete 
library (viridis) ## Cargar el paquete 

boxplot (expre, col = mako (6)) ## Boxplot de la base de datos con los perfiles de expresion en las distintas variantes o cepas


## Expresion diferencial 

types <- factor (c ("KO", "KO", "KO", "WT", "WT", "WT")) ## Dar la categoria de factor a al vector de...cepas u objetos de estudio
types ## Ver los "tipos", en este caso, solo con KO y WT, aunque, hay 6 en total, no toma los repetidos 
help ("factor")

help ("design")
design <- model.matrix (~ 0 + types) ## Crea una matriz a partir del objeto "types"
colnames (design) = levels(types) ## Los nombres de las columnas seran los factores de "types"
design ## Ver matriz

contMatrix = makeContrasts (KO - WT, levels = design) ## Usando la matriz anteriormente hecha, se hace un contraste o comparacion de
# KO con WT
contMatrix ## Ver la matriz de contraste


help ("lm.fit")
help ("contrasts.fit")
help ("eBayes")
fit  <- lmFit (expre, design) ## Modelo lineal de los datos originales y la matriz
fit2 <- contrasts.fit (fit, contMatrix) ## Modelo lineal del objeto "fit" (que ya es un modelo lineal) con la matriz de comparacion
fit2 <- eBayes (fit2) ## Analisis Bayesiano 

help ("topTable")
topTable (fit2, number = 10, sort.by = "p") ## Tabla de los genes mas importantes del modelo lineal, se puede elegir un parametro de cuantos genes
# se quieren ver en la tabla y se elige tambien que estadistico se quiere usar para esto, en este caso, fue p, ya que es el "estandar"


## Anotando los datos de expresion

library (mouse4302.db) ## Cargar la liberia (anteriormente instalada)
mouse4302 () ## Ver los mapeos en la libreria

fit2 ## Fir2 con TODA la informacion del experimento

head (fit2$genes) ## En esta parte no entendi muy bien por que me arroja un NULL, mas no error
probes <- fit2$genes$ID
probes ## Lo mismo 


## Asignar a objetos los nombres, simbolos y e identificadores de los genes
descriptions <- mget (probes, mouse4302GENENAME)
symbols <- mget (probes, mouse4302SYMBOL)
entrezids <- mget (probes, mouse4302ENTREZID)


## Sobrescribir esta informacion en el objeto fit2
fit2$genes$EntrezID <- unlist (entrezids)
fit2$genes$Symbol <- unlist (symbols)
fit2$genes$Description <- unlist (descriptions)


head (fit2$genes) ## Ver la primera parte de los genes en fit2, es DEMASIADO grande, incluso con head

help ("volcanoplot")
topTable (fit2, number = 10, sort.by = "p") ## Tabla de los genes mas importantes del modelo lineal, se puede elegir un parametro de cuantos genes
# se quieren ver en la tabla y se elige tambien que estadistico se quiere usar para esto, en este caso, fue p, ya que es el "estandar"
volcanoplot (fit2, highlight = 10, names = fit2$genes$Symbol) ## Crea una grafica de volcan con los 10 genes (usando el nombre corto) respecto al "log fold change"

deTable <- topTable (fit2, number = nrow (fit2), lfc = log2 (1.5), p.value = 0.05) ## Tabla para ver qué genes cumplen las condiciones de fold-change, este valor, lo podemos elegir nosotros
dim (deTable) ## Ver numero de genes que, en este caso, tienen un fold-change de 1.5 o más con un valor de p de 0.05

fullTable <- topTable (fit2, number = nrow (fit2), sort.by = "logFC") ## Practicamente lo mismo, pero con TODOS los genes
dim (fullTable)


## Ejercicio

# Cuantas sondas predicen como diferencialmente expresadas?
deTable <- topTable (fit2, number = nrow (fit2), lfc = log2 (1.5), p.value = 0.05) ## Tabla para ver qué genes cumplen las condiciones de fold-change, este valor, lo podemos elegir nosotros
dim (deTable) ## Ver numero de genes que, en este caso, tienen un fold-change de 1.5 o más con un valor de p de 0.05
## 219

# Cuantas decrementan y cuantas aumentan su expresion en el KO?
# Se que es un heat map lo que debemos hacer, pero no entiendo con que matriz se debe realizar 
heatmap (design)

# Cuantos genes unicos hay en estas listas? (?unique)

