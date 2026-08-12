library(ggplot2)
#library(dplyr)
library(rstatix)
library(tidyverse)

#Separar MCV del dataset completo
MCV_BDD_1025 <- subset(BDD_Art_curada_CCV_MCV_281025, Treatment == "MCV")

# PONDERACIÓN CON CRITERIOS ESTABLECIDOS - WPI specific weights  ---------------------------------

### Configuración de la ponderación:
##definición de pesos y patógenos específicos basado en el scorecard (Tabla S13)
PESOS_ESPECIFICOS <- c(
    "C. auris" = 0.5,
    "C. albicans" = 0.2,
    "C. albicans ATCC" = 0.2)

##Peso total restante para las demás cepas (1 - suma de pesos específicos)
PESO_TOTAL_ESPECIFICO <- sum(PESOS_ESPECIFICOS)
PESO_RESIDUAL_TOTAL <- 1 - PESO_TOTAL_ESPECIFICO
cat(paste("Peso total específico:", round(PESO_TOTAL_ESPECIFICO,4)))
cat(paste("\nPeso total residual para otros patógenos:", round(PESO_RESIDUAL_TOTAL, 4), "\n"))

##Calcular media de la inhibición y pivotar a formato ancho
datos_ancho <- MCV_BDD_1025 %>% 
    #Limpiar/seleccionar columnas clave
    select(Isolate, Pathogen, Inhibition) %>% 
    #Calcular medias 
    group_by(Isolate, Pathogen) %>% 
    summarise(Inhibicion_media = mean(Inhibition, na.rm = TRUE), 
              .groups = 'drop') %>% 
    #Pivotea 
    pivot_wider(names_from = Pathogen, values_from = Inhibicion_media)
print("Datos de inhibición media por patógeno (Formato ancho):")
print(head(datos_ancho))



##Z-scores --------------------------------------------------------
###ESTANDARIZACIÓN DE Z-SCORES

##Identificar todas las columnas de patógenos excepto Isolate
cols_patogenos_todos <- names(datos_ancho)[-1]
#Aplicar scale para calcular el z-score
datos_zscore <- datos_ancho %>% 
    mutate(across(all_of(cols_patogenos_todos), scale, .names = 
                      "Z_{.col}"))

# Configuración de pesos y Cálculo el Índice de Rendimiento Ponderado (IRP)
patogenos_especificos_nombres <- names(PESOS_ESPECIFICOS)
patogenos_restantes_nombres <- cols_patogenos_todos[! 
                                                        cols_patogenos_todos %in% 
                                                        patogenos_especificos_nombres]
n_restantes <- length(patogenos_restantes_nombres)

#Calculo del peso por patógeno
if (n_restantes == 0) {
    PESO_RESTANTE_INDIVIDUAL <- 0
} else {
    PESO_RESTANTE_INDIVIDUAL <- PESO_RESIDUAL_TOTAL / n_restantes    
}

cat(paste("Número de patógenos restantes:", 
          round(PESO_RESTANTE_INDIVIDUAL, 4), "\n"))

#Calculo del IRP
datos_ponderados <- datos_zscore  %>%
    rowwise()  %>%
    mutate(IRP_PARCIAL_ESPECIFICO = 0,
           IRP_PARCIAL_RESTI = 0)

#sUMA PONDERADA DE PATÓGENOS CON PESO ESPECIFICO
for (patogeno in patogenos_especificos_nombres) {
    col_z <- paste0("Z_", patogeno)
    peso <- PESOS_ESPECIFICOS[patogeno]
    #añadir contribución de cada patógeno específico
    datos_ponderados <- datos_ponderados %>%
        mutate(IRP_PARCIAL_ESPECIFICO
               =IRP_PARCIAL_ESPECIFICO + (.data[[col_z]] * peso))
}

#sUMA PONDERADA DE PATÓGENOS RESTANTES
if(n_restantes > 0) {
    #sumar los z score de patógenos restantes
    datos_ponderados <- datos_ponderados %>%
        mutate(Z_RESTO_SUMA = 
                   sum(c_across(paste0("Z_", 
                                       patogenos_restantes_nombres)))) %>%
        #Aplicar el peso residual individual
        mutate(IRP_PARCIAL_RESTO = Z_RESTO_SUMA * PESO_RESTANTE_INDIVIDUAL)
}


## WPI total ---------------------------------------------------------------
##calcular el IRP total
datos_ponderados <- datos_ponderados %>%
    mutate(IRP_FINAL = 
               IRP_PARCIAL_ESPECIFICO + 
               IRP_PARCIAL_RESTO) %>%
    ungroup() %>%
    
    #seleccionar columnas y ordenar 
    select(Isolate, IRP_FINAL) %>%
    arrange(desc(IRP_FINAL)) %>%
    rename(IRP = IRP_FINAL)

##Imprimir RAANKING final
cat("\nRanking Final por Índice de Rendimiento Ponderado (WPI):\n")
print(datos_ponderados)

#Visualización del ranking

grafico_IRP2 <- datos_ponderados %>%
    mutate(Isolate = factor(Isolate, levels = datos_ponderados$Isolate)) %>% # Ordenar las barras
    
    ggplot(aes(x = Isolate, y = IRP, fill = IRP)) +
    geom_bar(stat = "identity", color = "black") +
    
    geom_text(aes(label = round(IRP, 2)), vjust = -0.5, size = 3.5) +
    
    labs(
        title = paste0("Weighted Performance Index (WPI) - Specific Weights by pathogen"),
        subtitle = expression(paste(italic(" C. auris: 50%"),  italic(", C. albicans: 20%") )),               
        x = "Isolate",
        y = "WPI (Weighted Z score)",
        fill = "WPI"
    ) +
    theme_minimal(base_size = 12) +
    scale_fill_gradient(low = "#b2df8a", high = "#1f78b4") + # Paleta de colores más profesional
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# visuaización del gráfico
print(grafico_IRP2)
png(file="WPI-MCV_2_res350.png", res=350, width=3500, height=2000)
grafico_IRP2
dev.off()


write.csv(datos_ponderados, "WPI_auris50albicans20.csv")