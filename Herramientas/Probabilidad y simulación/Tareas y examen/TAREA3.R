library(MASS)
library(survival)
library(fitdistrplus)
library(foreign)

# TAREA 3

# 1. Distribución del mínimo de una muestra Uniforme(0,1)
min_values <- numeric(1000)
for (i in 1:1000) {
  uniform_sample <- runif(5000, 0, 1)
  min_values[i] <- min(uniform_sample)
}

# Ajustar distribución Beta
beta_fit <- fitdist(min_values, "beta")
summary(beta_fit)
plot(beta_fit)

# Parámetros estimados
beta_param1 <- beta_fit$estimate[1]
beta_param2 <- beta_fit$estimate[2]

# Bondad de ajuste
beta_gof <- gofstat(beta_fit, fitnames = "beta")
print(beta_gof)

# Test de Kolmogorov-Smirnov
ks.test(min_values, "pbeta", beta_param1, beta_param2)

# Valores teóricos y simulados
expected_value <- beta_param1 / (beta_param1 + beta_param2)
variance_theoretical <- (beta_param1 * beta_param2) / (((beta_param1 + beta_param2)^2) * (beta_param1 + beta_param2 + 1))
mean_simulated <- mean(min_values)
variance_simulated <- var(min_values)

cat("Media teórica:", expected_value, "\nVarianza teórica:", variance_theoretical, "\n")
cat("Media simulada:", mean_simulated, "\nVarianza simulada:", variance_simulated, "\n")

# Probabilidades P(X < 0.1)
prob_theoretical <- pbeta(0.1, beta_param1, beta_param2)
prob_simulated <- mean(min_values < 0.1)
cat("Probabilidad teórica P(X < 0.1):", prob_theoretical, "\nProbabilidad simulada P(X < 0.1):", prob_simulated, "\n")


# 2. Intervalo de confianza al 95% para una variable (Abdomen del archivo SPSS)
abdomen_data <- read.spss("_Pallarés_Díez.sav", to.data.frame = TRUE)$Abdomen

# Ajustar distribución Normal
normal_fit <- fitdist(abdomen_data, "norm")
summary(normal_fit)
plot(normal_fit)

# Estadísticas de bondad de ajuste
normal_gof <- gofstat(normal_fit, fitnames = "norm")
print(normal_gof)

# Test de Kolmogorov-Smirnov
ks.test(abdomen_data, "pnorm", mean(abdomen_data), sd(abdomen_data))

# Calcular intervalo de confianza
mean_abdomen <- numeric(1000)
for (i in 1:1000) {
  mean_abdomen[i] <- mean(sample(abdomen_data, size = 100, replace = TRUE))
}

CI_abdomen <- quantile(mean_abdomen, c(0.025, 0.975))
cat("Intervalo de confianza al 95% para Abdomen:", CI_abdomen, "\n")


# 3. Probabilidad con Gamma(5, 2)
alpha_gamma <- 5
beta_gamma <- 2

# Cálculos teóricos con TCL
mean_gamma <- alpha_gamma / beta_gamma
variance_gamma <- alpha_gamma / (beta_gamma^2) / 1000
probability_tcl <- 1 - pnorm(3, mean = mean_gamma, sd = sqrt(variance_gamma))
cat("Probabilidad teórica con TCL P(X̄ > 3):", probability_tcl, "\n")

# Simulaciones con Gamma
mean_values_gamma <- numeric(1000)
for (i in 1:1000) {
  gamma_sample <- rgamma(1000, shape = alpha_gamma, rate = beta_gamma)
  mean_values_gamma[i] <- mean(gamma_sample)
}

probability_simulated <- mean(mean_values_gamma > 3)
cat("Probabilidad simulada P(X̄ > 3):", probability_simulated, "\n")

