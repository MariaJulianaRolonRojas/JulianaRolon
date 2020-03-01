#CONJUNTO FX PARA DESVEST Y NORMALIZACIÓN DE LECTURA 
#install.packages("ggplot2")
library("ggplot2")
mediaydesvest <- function(filename,columnas,blanco_multi=FALSE, lineas=1) { #Numero de columnas que el vittor analiza
  print(filename)
  #Transforma archivo en vector de lineas 
  lines=scan(filename, what="character",sep="\n")
  print(lines[1])
  #Grep Busca ese Texto dentro de la lista
  divisiones=grep("Plate	Repeat	End time	Start temp.	End temp.	BarCode",lines)
  print(divisiones)
  #Transforma listas de la tabla en documento de texto y lo lee como tabla hasta donde yo le diga
  dat=read.table(textConnection(lines[1:(divisiones[1]-1)]), header=TRUE,sep="\t")
  #Encontrar ultimo poso 
  x=tail(dat$Well,n=1)
  #Coge la letra del poso y la vuelve numero para saber cuantos posos hay en el experimento
  letra=utf8ToInt(substr(x, start = 1, stop = 1))-utf8ToInt("A")+1
  #numero despues de la letra en la lectura para numero de posos
  numero=strtoi(substr(x, start = 2, stop = 3))
  #calcula el numero de posos utilizados
  posos=(letra-1)*columnas+numero
  #Encuentra el numero de lecturas que se hicieron en el experimento
  repeticion=tail(dat$Repeat,n=1)
  
  #Encuentra lo que calcula Abs, GFP, Cherry...
  dims = (length(dat)-4)/2
  #Crea matriz para almacenar los datos 
  result=array(dim=c(dims,repeticion,posos))
  resultName=c()
  #Pasa por todas las dimensiones de la tabla y los llena
  for(j in 1:dims){
    pos = 4+j*2  #posicion de la columna
    #pasa por todas las repeticiones y pone los datos
    for (i in 1:repeticion) {	
      result[j,i,]=dat[((i-1)*posos+1):(i*posos),pos]
    }
    #Guarda el nombre de la columna como Abs,GFP,Cherry
    colName = strsplit(names(dat)[pos],"\\.")[[1]][1]
    resultName = c(resultName, colName)
  }
  #print(dim(result))
  #print(result[1,,])
  
  #matriz con los datos y les resta el blanco de cada lectura
  nblanks=1
  if(blanco_multi){
    nblanks=4
  }
  background=array(dim=c(dims,repeticion,posos-nblanks*3))
  if(blanco_multi){
    #blanco multiple
    nonblank=posos-(3*nblanks)
    blank=(nonblank/3)+1
    print(blank)
    for(j in 1:dims){
      #pasa por todas las repeticiones y pone los datos 
      for (i in 1:repeticion) {
        k=blank
        for(p in 1:nonblank){
          background[j,i,p]=result[j,i,p]-mean(result[j,i,((k-1)*3+1):(k*3)])
          if(p %% 3 == 0){
            if(k<blank+3){
              k=k+1
            }else{
              k=blank
            }
          }
        }
      }
    }
  }
  else{
    #blanco unico
    for(j in 1:dims){
      #pasa por todas las repeticiones y pone los datos 
      for (i in 1:repeticion) {
        background[j,i,1:(posos-3)]=result[j,i,1:(posos-3)]-mean(result[j,i,(posos-2):posos])
      }
    }
  }
  #print(background[1,,])
  
  #Fluorecencia Normalizada
  ByOD=array(dim=c(dims-1,repeticion,posos-3*nblanks))
  for(j in 2:dims){
    ByOD[j-1,,]=background[j,,]/background[1,,]
  }
  #print(ByOD[1,,])
  
  #Calcula las triplicatas
  ntri = (posos-3*nblanks)/(3*lineas)
  media=array(dim=c(dims,repeticion,ntri))
  desvest=array(dim=c(dims,repeticion,ntri))
  #Loop para abs,GFP,mCherry
  for(k in 1:dims){
    # #Loop para cada lectura
    # for (j in 1:repeticion){
    #   #Loop para cada triplicata
    #   for (i in 1:((posos-3*nblanks)/3)){
    #     if(k==1){
    #       #guarda la media de las triplicatas del od 
    #       dados = background[1,j,((i-1)*3+1):(i*3)]
    #     } else {
    #       #guada la media de las triplicatas de las flourescencias normalizadas
    #       dados = ByOD[(k-1),j,((i-1)*3+1):(i*3)]
    #     }
    #     media[k,j,i] = mean(dados)
    #     desvest[k,j,i] = sd(dados)
    #   }
    # }
    #Loop para cada lectura
    for (j in 1:repeticion){
      nl = 1
      for(i in 1:ntri){
        ind = (((i-1)*3+1) %% 12)+(nl-1)*12
        dados = c()
        for(l in 1:lineas){
          interv = (ind+(l-1)*12):(ind+2+(l-1)*12)
          if(k==1){
            #guarda la media de las triplicatas del od 
            dados = c(dados, background[1,j,interv])
          } else {
            #guada la media de las triplicatas de las flourescencias normalizadas
            dados = c(dados, ByOD[(k-1),j,interv])
          }
        }
        media[k,j,i] = mean(dados)
        desvest[k,j,i] = sd(dados)
        if((ind+2) %% 12==0){
          nl=nl+lineas
        }
      }
    }
  }
  return(list(media=media,desvest=desvest,tipos=resultName))
}

toggplot <- function(triplicata, concentraciones, compuestos, lectura=0.5){
  #Los datos como vienen estan en 3 dimensiones, hay que combiarlos a tabla entonces
  #Crea una tabla con los datos de dimensiones de las medias y el df2 para desvest
  df1 = as.data.frame.table(triplicata$media, responseName = "media")
  df2 = as.data.frame.table(triplicata$desvest, responseName = "sd")
  #Lo anterior creo una tabla para cada, este comando lo une. Como? crea una columna y le atribuye los datos
  df1$sd = df2$sd
  #Vuelve numericos los valores de las repeticiones y de las muestras
  df1$Var2 = as.numeric(as.factor(df1$Var2))*lectura
  df1$Var3 = as.numeric(as.factor(df1$Var3))
  #Se crean las concentraciones y el compuesto CAMBIA CAMBIA CAMBIA AAAAAQUIIII <--------------------------
  #Crea columna vacia y le pone el nombre
  df1$conc = NA
  df1$comp = NA
  df1$ctrl = NA
  a = 1
  b = 1
  ctrl=FALSE
  for (i in seq(1,dim(triplicata$media)[3])){
    df1[df1$Var3==i,]$conc = concentraciones[a]
    df1[df1$Var3==i,]$comp = compuestos[b]
    df1[df1$Var3==i,]$ctrl = ctrl
    a=a+1
    if(a>length(concentraciones)){
      b=b+1
      a=1
    }
    if(b>length(compuestos)){
      b=1
      ctrl=TRUE
    }
  }
  #Cambia las letras extrañas que salian por el nombre del tipo
  df1$Var1 = as.numeric(as.factor(df1$Var1))
  for (i in 1:length(triplicata$tipos)){
    df1$Var1 = replace(df1$Var1, df1$Var1==i, triplicata$tipos[i])
  }
  return(df1)
}

plotConcentracion = function(plot, tipo, compound, ctrl, lconc=c("All")){
  if(lconc==c("All")){
    plotOD = plots[plots$Var1==tipo & plots$comp==compound & plots$ctrl==ctrl,]
  }else{
    plotOD = plots[plots$Var1==tipo & plots$comp==compound & plots$ctrl==ctrl & plots$conc %in% lconc,]
  }
  
  if(tipo=="Absorbance"){
    titulo=paste("Bacterial growth in",compound)
    ytext="OD"
  }else{
    titulo=paste("Bacterial expression of", tipo, "in", compound)
    ytext=paste(tipo,"/OD")
  }
  
  if(ctrl){
    titulo=paste(titulo,"(control)")
  }
  
  p = ggplot(plotOD, aes(x=Var2, y=media, group=conc, color=conc)) + 
    geom_errorbar(aes(ymin=media-sd, ymax=media+sd), width=.1) +
    geom_line() + geom_point()+
    scale_color_brewer(palette="Spectral")+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p+labs(title=titulo, x="Hours", y=ytext, color="Concentrations")+ # for the main title, axis labels and legend titles
    theme(plot.title = element_text(hjust = 0.5))+
    theme(legend.position=c(0.2,0.8))
}

#Busca Archivo, AQUI COMIENZA EL CODIGO! CORRERLO SOLO DESDE AQUI 
filename1=file.choose()
#calcula media, desvio, quita el ruido de fondo del blanco de todos los posos en todas las lecturas 
triplicata = mediaydesvest(filename1,12,blanco_multi=TRUE,lineas=3)

#Guarda los datos en la carpeta que le ponga
folder = ""
for(k in 1:dims){
  dirout = paste0(folder,resultName[k],".txt")
  write.table(triplicata[k,,],file=dirout, sep="\t",row.names=FALSE)
}


concentraciones = c("0mM","0.001mM","0.01mM","0.1mM")
#compuestos = c("Ferulic Acid","Vanillic Acid", "Vanillin")
compuestos = c("Vanillin")

#Transforma los datos de 3 dimensiones para 2 dimensiones para crear una tabla y sea mas facil hacer el grafico
plots = toggplot(triplicata, concentraciones, compuestos)

#Grafico

#Plot de Expresion/Crecimiento Vs Lectura en Concentraciones de compuestos (este me dice cuanto expreso para cada compuesto y si se inhibio el crecimiento con laguno de ellos)
plotConcentracion(plots, tipo = "mCherry", compound = "Vanillin", ctrl = FALSE)
plotConcentracion(plots, tipo = "mCherry", compound = "Vanillin", ctrl = FALSE, lconc = c("0mM","0.01mM","0.1mM"))

for(ctrl in c(FALSE,TRUE)){
  for(t in 1:length(triplicata$tipos)){
    for (c in 1:length(compuestos)) {
      print(plotConcentracion(plots, tipo =triplicata$tipos[t], compound = compuestos[c], ctrl))
    } 
  }
}

#plots sem 0.001mM
for(ctrl in c(FALSE,TRUE)){
  for(t in 1:length(triplicata$tipos)){
    for (c in 1:length(compuestos)) {
      print(plotConcentracion(plots, tipo =triplicata$tipos[t], compound = compuestos[c], ctrl, lconc = c("0mM","0.01mM","0.1mM")))
    } 
  }
}





