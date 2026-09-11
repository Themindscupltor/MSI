#!/usr/bin/Rscript
##MarcoG
#Analisis demografico - GBE
#Ultima fecha de edición: 10/sep/2026

#Establecemos la seed
set.seed(seed = 123456789)

# Cargando las librerias
library(dartR.popgen)
library(dartRverse)
library(vcfR)

#Establecemos para empezar el inicio de la corrida del codigo
fechaInicio<-date()

# Establecemos el directorio de trabajo
#setwd(dir="/run/media/pinkpunk/myriad/GBE/starway") <- directorio original
setwd(dir = "/run/media/pinkpunk/myriad/GBE/starway")

#------------------------------------
# Carga de datos
#------------------------------------

#cargamos la informacion del hapmap
hapmap<-read.table(file = "./hapmap.tsv",
                   header = TRUE,
                   sep = "\t")

# Cargamos el archivo VCF completo
asty <- gl.read.vcf(
  vcffile="./Amex_rad-snpsFiltered_m01_ld06.vcf.gz",
  verbose=5
)

#---------------------------------------------------
# Quality check 0w0
#---------------------------------------------------

# Checamos la compatibilidad con dartR
asty <- gl.compliance.check(
  x=asty,
  verbose=5
)

#Agregamos la informacion por poblacion
pop(asty)<-hapmap$poblacion

#Hacemos la imputacion basado enla frecuencia
#alelica poblacional
asty_imp <- gl.impute(
  x = asty,
  method = "frequency",
  verbose = 5
)

# ---------------------------------------------------------
# ANALISIS DEMOGRAFICO  - INICIO 
# ---------------------------------------------------------
#Creamos una lista donde guardaremos nuestros
#Resultados
resultados<-list()

# Stairway2 genera el script de ejecucion y sus archivos
# auxiliares en tempdir(). El script utiliza rutas relativas,
# por lo que debe ejecutarse desde ese directorio.

#Vemos donde esta nuestro directorio temporal
#En este se correran los analisis demograficos
tempdir()

#NOTA: aca vamos a homogenizar las pruebas
#Normalmente el analisis se hace sobre
#snps distintos en cada corrida, hasta ahora
#vamos a generar un vector con las posiciones en la matriz
#para que cada corrida se haga sobre los mismos sitios
sitios<-sample(x = 1:139097,size = 1000,replace = FALSE) %>%
  sort()
#Eliminamos los contenidos de la carpeta
#pues el contenido puede afectar la corrida
#del programa de stairway
unlink(list.files(path = tempdir(), 
                  full.names = TRUE), 
       recursive = TRUE)

#Revisamos que no haya nada aca.
list.files(tempdir())

#--------------------------------------
#  Stairway plot2 - Cuevas L1
#--------------------------------------
#para poblaciones de cueva de L1
#sin considerar Subterraneo ni los individuos de Arroyo
filtro <- hapmap %>%
  subset(lineage=="L1" & 
         environment == "Cave" &
          poblacion!=c("Subterraneo") &
           poblacion!=c("Arroyo"))

#Revisamos 
#Atento a que coincidan estos individuos
asty_filtrado<-asty_imp[asty_imp@ind.names  %in% filtro$id,]
asty_filtrado@ind.names

# Descargamos Stairway2 en el directorio temporal de R
gl.download.binary(
  software="stairway2",
  os="linux"
)

#Hacemos el modelado demografico sobre las poblaciones
#En este primer caso, cuevas de Linaje 1
caveL1<-gl.run.stairway2(x = asty_filtrado[,sitios],
                         L = 129100,
                         mu=3.5e-9,
                         gentime=1,
                         blueprint="caveL1_starway",
                         filename="caveL1_starway",
                         pct_training=0.70,
                         seed=12345,
                         stairway2.path=file.path(tempdir(),"stairway2"),
                         parallel=4,
                         nreps=200,
                         run=TRUE,
                         cleanup = FALSE
)

#Guardamos los resultados 
resultados$CaveL1 <- caveL1$history %>%
  mutate(label = "Cave L1")

#Eliminamos los contenidos de la carpeta
unlink(list.files(path = tempdir(), 
                  full.names = TRUE), 
       recursive = TRUE)

#--------------------------------------
#  Stairway plot2 - Cuevas L2
#--------------------------------------
#Volvemos a descargar el programa
gl.download.binary(
  software="stairway2",
  os="linux"
)

#generemos los resultados de la segunda corrida
filtro <- hapmap %>%
  subset(lineage=="L2" & 
           environment == "Cave" &
           poblacion != "Arroyo" &
           poblacion != "Chica" &
           poblacion != "Toro" )

#Revisamos 
#Atento a que coincidan estos individuos
asty_filtrado<-asty_imp[asty_imp@ind.names  %in% filtro$id,]
asty_filtrado@ind.names

#Hacemos el modelado demografico
caveL2<-gl.run.stairway2(x = asty_filtrado[,sitios],
                         L = 129100,
                         mu=3.5e-9,
                         gentime=1,
                         blueprint="caveL2_starway",
                         filename="caveL2_starway",
                         pct_training=0.70,
                         seed=12345,
                         stairway2.path=file.path(tempdir(),"stairway2"),
                         parallel=4,
                         nreps=200,
                         run=TRUE,
                         cleanup = FALSE
)


#Agregamos los resultados de la segunda corrida
resultados$CaveL2 <- caveL2$history %>%
  mutate(label = "Cave L2")

#Eliminamos los contenidos de la carpeta
unlink(list.files(path = tempdir(), 
                  full.names = TRUE), 
       recursive = TRUE)

#--------------------------------------
#  Stairway plot2 - Superficies L1
#--------------------------------------
#Volvemos a descargar el programa
gl.download.binary(
  software="stairway2",
  os="linux"
)

#generemos los resultados de la segunda corrida
filtro <- hapmap %>%
  subset(lineage=="L1" & 
           environment == "Surface" &
           poblacion!="Caballo Moro_S")
#Revisamos 
#Atento a que coincidan estos individuos
asty_filtrado<-asty_imp[asty_imp@ind.names  %in% filtro$id,]
asty_filtrado@ind.names

#Hacemos el modelado demografico sobre las poblaciones
#En este primer caso, cuevas de Linaje 1
surfL1<-gl.run.stairway2(x = asty_filtrado[,sitios],
                         L = 129100,
                         mu=3.5e-9,
                         gentime=1,
                         blueprint="surfL1_starway",
                         filename="surfL1_starway",
                         pct_training=0.70,
                         seed=12345,
                         stairway2.path=file.path(tempdir(),"stairway2"),
                         parallel=4,
                         nreps=200,
                         run=TRUE,
                         cleanup = FALSE
)

#Agregamos los resultados de la segunda corrida
resultados$surfL1 <- surfL1$history %>%
  mutate(label = "Surface L1")

#Eliminamos los contenidos de la carpeta
#pues pueden estar alterando las llamadas
#del programa de stairway
unlink(list.files(path = tempdir(), 
                  full.names = TRUE), 
       recursive = TRUE)

#--------------------------------------
#  Stairway plot2 - Superficies L2
#--------------------------------------
#Volvemos a descargar el programa
gl.download.binary(
  software="stairway2",
  os="linux"
)

#generemos los resultados de la segunda corrida
filtro <- hapmap %>%
  subset(lineage=="L2" & 
           environment == "Surface")

#Revisamos 
#Atento a que coincidan estos individuos
asty_filtrado<-asty_imp[asty_imp@ind.names  %in% filtro$id,]
asty_filtrado@ind.names

#Hacemos el modelado demografico sobre las poblaciones
surfL2<-gl.run.stairway2(x = asty_filtrado[,sitios],
                         L = 129100,
                         mu=3.5e-9,
                         gentime=1,
                         blueprint="surfL2_starway",
                         filename="surfL2_starway",
                         pct_training=0.70,
                         seed=12345,
                         stairway2.path=file.path(tempdir(),"stairway2"),
                         parallel=4,
                         nreps=200,
                         run=TRUE,
                         cleanup = FALSE
)

#Agregamos los resultados de la segunda corrida
resultados$surfL2 <- surfL2$history %>%
  mutate(label = "Surface L2")

#Eliminamos los contenidos de la carpeta
#pues pueden estar alterando las llamadas
#del programa de stairway
unlink(list.files(path = tempdir(), 
                  full.names = TRUE), 
       recursive = TRUE)

#--------------------------------------------------------
# Stairway plot2 - Apartado grafico manual
#--------------------------------------------------------
#vamos a hacer un rbind entre los dos
#sets de datos conjuntos

#Juntamos todos los datos
datosConjuntos <- rbind.data.frame(resultados$surfL2,
                                   resultados$surfL1,
                                   resultados$CaveL1,
                                   resultados$CaveL2
                                   )

#Graficamos utilizando ggplot2
ggplot(data = datosConjuntos, 
       mapping = aes(x = year, y = Ne_median)) + 
  geom_point(mapping = aes(color=label)) +
  facet_wrap(~label)+
  geom_line() + 
  geom_ribbon(aes(ymin = low95, 
                  ymax = high95), 
              alpha = 0.2) +
  ylab("Effective population size") + 
  xlab("Time (1 year per gen.)")+
  xlim(NA,25000) + 
  ylim(NA,4e5) +
  theme_bw() + 
  theme(legend.position = "none")

#Vamos a guardar los resultados generados 
write.csv(x = datosConjuntos,
          file = "./resultados_STplot-test.csv",
          sep = ",")

#Imprimimos la estampa de la fecha
fechaFin<-date()
cat("inicio: ",fechaInicio,"\n","fin: ",fechaFin, "\n")

#------------------------
# Fin de codigo
#------------------------
#    Marco G  0w0 ⠀⠀⠀⠀⠀⠀⠀⠀⣠⣄⣠⣄⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⠏⠀
# ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⣀⣀⡀⠀⠀⠀⠀⠀⠀⠈⠀⠀⠀
# ⠀⠀⠀⠀⣀⡤⣤⠶⠛⠉⠉⠀⠀⠉⠉⠛⠲⣤⣤⣄⠀⠀⠀⠀⠀
# ⠀⠀⠀⡼⠃⠈⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠁⠈⢧⠀⠀⠀⠀
# ⠀⠀⡼⢁⡆⠀⢀⡀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡄⠀⢰⡈⢧⠀⠀⠀
# ⢀⡞⠁⣸⠁⠀⢠⠬⠓⠀⣀⠀⠀⢀⠀⠛⠥⣄⠀⠀⡇⠈⢳⡀⠀
# ⡞⠀⠀⠹⣆⢰⣒⠆⠀⠀⠓⠊⠙⠚⠁⠀⠸⠭⠇⣰⠇⠀⠀⢻⠀
# ⣇⠀⠀⠀⠈⣹⠶⠦⠤⠤⣤⣤⣤⡤⠤⠤⠴⠶⣏⠁⠀⠀⠀⣸⠁
# ⠈⠓⠲⠖⠚⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠓⠲⠶⠚⠁⠀
