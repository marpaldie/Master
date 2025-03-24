install.packages('fitdistrplus')
install.packages("survival")
install.packages("readr")
install.packages("ggplot2")
install.packages("univariateML")
install.packages("tidyverse")

rm(list = ls())

# Cargar las librerías necesarias
library(tidyverse)
library(readr)
library(MASS)
library(survival)
library(fitdistrplus)
library(ggplot2)
library(univariateML)
library(dplyr)

# Cargar la librería necesaria
library(fitdistrplus)

# Leer los datos
Super <- read.csv("C:/Users/Maria/Desktop/UV/Herramientas_matemáticas_informáticas/PROBABILIDAD/R/Super.csv", sep="")
View(Super)
head(Super)
Temps <- Super$Temps

# Ajustar las distribuciones
fgamma <- fitdist(Temps, "gamma")           
fexp <- fitdist(Temps, "exp")               
fnorm <- fitdist(Temps, "norm")             
fchi <- fitdist(Temps, "chisq", start=list(df=2)) 
funif <- fitdist(Temps, "unif")             
fweibull <- fitdist(Temps, "weibull")  
flognorm <- fitdist(Temps, "lnorm")  

# Crear el gráfico de comparación
par(mfrow=c(2,2))
plot.legend <- c("Gamma", "Exponential", "Normal", "Chi-Squared", "Uniform", "Weibull", "Log-normal")

# Comparación de distribuciones usando los gráficos de `fitdistrplus`
denscomp(list(fgamma, fexp, fnorm, fchi, funif, fweibull, flognorm), legendtext=plot.legend)  # Comparación de densidades
qqcomp(list(fgamma, fexp, fnorm, fchi, funif, fweibull, flognorm), legendtext=plot.legend)    # Comparación QQ plot
cdfcomp(list(fgamma, fexp, fnorm, fchi, funif, fweibull, flognorm), legendtext=plot.legend)   # Comparación CDF
ppcomp(list(fgamma, fexp, fnorm, fchi, funif, fweibull, flognorm), legendtext=plot.legend)    # Comparación PP plot

summary(fgamma)
summary(fexp)
summary(fnorm)
summary(fchi)
summary(funif)
summary(fweibull)
summary(flognorm)

#Representamos log-normal que es la que tiene mejores valors AIC y BIC
par(mfrow=c(2,2))
plot.legend <- c("Log-normal")
denscomp(list(flognorm), legendtext=plot.legend)  # Comparación de densidades
qqcomp(list(flognorm), legendtext=plot.legend)    # Comparación QQ plot
cdfcomp(list(flognorm), legendtext=plot.legend)   # Comparación CDF
ppcomp(list(flognorm), legendtext=plot.legend)    # Comparación PP plot

# Calcular la probabilidad de supervivencia para 12 y 24 meses con la distribución Log-normal
prob_12 <- 1 - plnorm(12, meanlog = flognorm$estimate[1], sdlog = flognorm$estimate[2])
prob_24 <- 1 - plnorm(24, meanlog = flognorm$estimate[1], sdlog = flognorm$estimate[2])

# Mostrar las probabilidades calculadas
cat("Probabilidad de supervivencia de al menos 12 meses:", prob_12, "\n")
cat("Probabilidad de supervivencia de al menos 24 meses:", prob_24, "\n")

# Crear una secuencia de tiempos para graficar
time_values <- seq(0, max(Temps), length.out = 100)

# Calcular la función de supervivencia para cada tiempo
surv_values <- 1 - plnorm(time_values, 
                          meanlog = flognorm$estimate["meanlog"], 
                          sdlog = flognorm$estimate["sdlog"])

# Crear un dataframe para la gráfica
surv_data <- data.frame(Time = time_values, Survival_Probability = surv_values)

# Graficar la función de supervivencia
ggplot(surv_data, aes(x = Time, y = Survival_Probability)) +
  geom_line(color = "red", size = 1) +
  labs(title = "Función de Supervivencia - Log-normal",
       x = "Tiempo (meses)",
       y = "Probabilidad de Supervivencia") +
  scale_x_continuous(breaks = seq(0, max(time_values), by = 15)) +  # Marcas del eje X cada 15 meses
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +  # Marcas del eje Y cada 0.1 (10%)
  theme_minimal()

# Crear una secuencia de tiempos para graficar la densidad de probabilidad
time_values <- seq(0, max(Temps), length.out = 100)

# Calcular la función de densidad de probabilidad para cada tiempo
meanlog <- flognorm$estimate["meanlog"]  
sdlog <- flognorm$estimate["sdlog"]      

# Calcular la función de densidad de probabilidad para cada tiempo usando los parámetros ajustados por `fitdist`
density_values <- dlnorm(time_values, meanlog = meanlog, sdlog = sdlog)

# Crear un dataframe para la gráfica
density_data <- data.frame(Time = time_values, Density = density_values)

# Graficar la función de densidad de probabilidad
ggplot(density_data, aes(x = Time, y = Density)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "Función de Densidad de Probabilidad - Log-normal",
       x = "Tiempo (meses)",
       y = "Densidad de Probabilidad") +
  theme_minimal()



