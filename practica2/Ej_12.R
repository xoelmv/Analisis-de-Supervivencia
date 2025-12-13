# Ejercicio 12.
# 
# Considérense los n = 1835 adultos (edad ≥ 18) con leucemia
# en el conjunto de datos ebmt1{mstate}, que pueden experimentar una recaída
# (evento 1) o una muerte en remisión (evento 2) a lo largo del seguimiento.
# Se pide:
# 
# 12.1 Mediante un modelo de regresión de riesgos de subdistribución proporcionales,
#      estimar el efecto del score. ¿Es estadísticamente significativo al 5%?
#      Obtener los factores multiplicativos de riesgo para los riesgos de subdistribución
#      tomando el score de riesgo bajo como grupo de referencia.
# 
# 12.2 Discutir las limitaciones o posibles inconsistencias del modelo utilizado.
# 
# 12.3 Representar gráficamente las funciones de incidencia acumulada estimadas
#      a partir del modelo para los tres grupos de score.
# 
# 12.4 Construir un gráfico de bondad de ajuste para el modelo utilizado.

library(mstate)

data(ebmt1)

# ------------------------------------------------------------------------------
# Preparamos los datos
# ------------------------------------------------------------------------------

ebmt1_adultos <- ebmt1[ebmt1$age>= 18,] # filtramos para quedarnos con los adultos
# dim(ebmt1_adultos)
# head(ebmt1_adultos, 10)
# summary(ebmt1_adultos)

# si rel < srv -> recaida
ebmt1_adultos$ftime <- pmin(ebmt1_adultos$rel, ebmt1_adultos$srv)

# variable de estado:
# 0 = censura
# 1 = recaída (relstat == 1 y rel < srv)
# 2 = muerte en remisión (relstat == 1 y rel == srv)
ebmt1_adultos$fstatus <- ifelse(ebmt1_adultos$relstat == 1 & ebmt1_adultos$rel < ebmt1_adultos$srv, 1,
                         ifelse(ebmt1_adultos$srvstat == 1 & ebmt1_adultos$rel == ebmt1_adultos$srv, 2, 
                                    0)) 

table(ebmt1_adultos$fstatus)
# 773 observaciones censuradas
# 421 sufieron el evento 1 (Recaída)
# 641 sufrieron el evento 2 (Muerte en remisión)

datos <- data.frame(id=ebmt1_adultos$patid,
                    time=ebmt1_adultos$ftime,
                    status=ebmt1_adultos$fstatus,
                    score = ebmt1_adultos$score,
                    age = ebmt1_adultos$age)

# Crear variables dummy para el score
scoreMedium <- ifelse(datos$score == "Medium risk", 1, 0)
scoreHigh <- ifelse(datos$score == "High risk", 1, 0)

score <- cbind(scoreMedium, scoreHigh)
head(score)



# ------------------------------------------------------------------------------
# Ajustamos el modelo de regresión de riesgos de subdistribución proporcionales
# ------------------------------------------------------------------------------
library(cmprsk) 
modelo_1 <- crr(ftime = datos$time,
              fstatus = datos$status,
              cov1 = score,
              failcode=1, # recaida
              cencode=0)

summary(modelo_1)

# La variable `score` es significativa. Los coeficientes de ambos niveles del score
# presentan p-valores menores a 0.05:
#   -  scoreMedium: p-value = 0.0440
#   -  scoreHigh: p-value = 0.0018
# 
# Factores multiplicativos de riesgo para los riesgos de subdistribución tomando
# el score de riesgo bajo como grupo de referencia
#   - scoreMedium: exp(0.277) = 1.32
#   - scoreHigh: exp(0.604) = 1.83


# Se ajustó un modelo de riesgos de subdistribución proporcionales para evaluar el
# efecto del score sobre el riesgo de recaída (evento 1), considerando la muerte en
# remisión como evento competitivo. 

# Los resultados muestran que el efecto del score es estadísticamente significativo 
# al 5\% (p = 0.044 para "Medium risk" y p = 0.0018 para "High risk").
# Tomando como referencia el grupo de bajo riesgo ("Low Risk"), el grupo de riesgo medio
# tiene un 32% mayor riesgo de recaída y el grupo de alto riesgo un 83\% mayor riesgo.

# Se ajustó un segundo modelo para estudiar el efecto del score sobre la muerte en remisión,
# considerando la recaída como evento competitivo.

modelo_2 <- crr(ftime = datos$time,
                fstatus = datos$status,
                cov1 = score,
                failcode=2, # muerte en remisión
                cencode=0)

summary(modelo_2)

# El grupo de riesgo medio tiene 1.86 veces más riesgo subdistribucional de morir en remisión.

# El grupo de alto riesgo casi triplica el riesgo (HR = 2.86) comparado con el grupo de bajo riesgo.

# Ambos resultados son altamente significativos (p < 0.001).


# 12.2
# En modelos de Fine & Gray, cada evento se modela por separado, especificando cuál
# es el evento de interés y tratando los demás como riesgos competitivos. Por tanto,
# los resultados de cada modelo solo pueden interpretarse dentro de su contexto 
# específico (recaída o muerte).

# -------------------------------------------------------------------------
# Funciónes de Incidencia Acumulada (Cumulative Incidence Function, CIF)
# -------------------------------------------------------------------------

par(mfrow=c(1,2))

# CIF para "Low" (categoría de referencia: scoreMedium = 0, scoreHigh = 0)
pred_low_1 <- predict(modelo_1, cov1 = matrix(c(0, 0), nrow = 1))
plot(pred_low_1, ylim = c(0, 1), col="green", xlab = "Tiempo", ylab = "Probabilidad de sufrir recaída", main = "CIFs según Score")

# CIF para "Medium" (scoreMedium = 1, scoreHigh = 0)
pred_medium_1 <- predict(modelo_1, cov1 = matrix(c(1, 0), nrow = 1))
lines(pred_medium_1, lty=2, col="blue")

# CIF para "High" (scoreMedium = 0, scoreHigh = 1)
pred_high_1 <- predict(modelo_1, cov1 = matrix(c(0, 1), nrow = 1))
lines(pred_high_1, lty=3, col="red")

# Leyenda
legend("topright", legend = c("Low", "Medium", "High"),
       lty = 1:3, col=c("green", "blue", "red"))


# CIF para "Low" (categoría de referencia: scoreMedium = 0, scoreHigh = 0)
pred_low_2 <- predict(modelo_2, cov1 = matrix(c(0, 0), nrow = 1))
plot(pred_low_2, ylim = c(0, 1), col="green", xlab = "Tiempo", ylab = "Probabilidad de muerte en remisión", main = "CIFs según Score")

# CIF para "Medium" (scoreMedium = 1, scoreHigh = 0)
pred_medium_2 <- predict(modelo_2, cov1 = matrix(c(1, 0), nrow = 1))
lines(pred_medium_2, lty=2, col="blue")

# CIF para "High" (scoreMedium = 0, scoreHigh = 1)
pred_high_2 <- predict(modelo_2, cov1 = matrix(c(0, 1), nrow = 1))
lines(pred_high_2, lty=3, col="red")

# Leyenda
legend("topright", legend = c("Low", "Medium", "High"),
       lty = 1:3, col=c("green", "blue", "red"))




# ------------------------------------------------------------------------------
# 12.1.1 Modelo Fine-Gray para recaída (evento 1)
# ------------------------------------------------------------------------------

# Ajustamos el modelo de riesgos de subdistribución para recaída
fg_recaida <- crr(ftime = datos$time,
                  fstatus = datos$status,
                  cov1 = model.matrix(~ score, data = datos)[,-1], # Eliminamos intercept
                  failcode = 1,
                  cencode = 0)

# Resumen del modelo
summary_fg_recaida <- summary(fg_recaida)
print(summary_fg_recaida)

# Coeficientes y significación estadística
cat("\nCoeficientes para recaída (evento 1):\n")
print(summary_fg_recaida$coef)

# Factores multiplicativos de riesgo (HR)
cat("\nFactores multiplicativos de riesgo (HR) para recaída:\n")
cat("Medium vs Low risk:", exp(summary_fg_recaida$coef[1,1]), "\n")
cat("High vs Low risk:", exp(summary_fg_recaida$coef[2,1]), "\n")

# ------------------------------------------------------------------------------
# 12.1.2 Modelo Fine-Gray para muerte en remisión (evento 2)
# ------------------------------------------------------------------------------

# Ajustamos el modelo de riesgos de subdistribución para muerte
fg_muerte <- crr(ftime = datos$time,
                 fstatus = datos$status,
                 cov1 = model.matrix(~ score, data = datos)[,-1],
                 failcode = 2,
                 cencode = 0)

# Resumen del modelo
summary_fg_muerte <- summary(fg_muerte)
print(summary_fg_muerte)

# Coeficientes y significación estadística
cat("\nCoeficientes para muerte en remisión (evento 2):\n")
print(summary_fg_muerte$coef)

# Factores multiplicativos de riesgo (HR)
cat("\nFactores multiplicativos de riesgo (HR) para muerte:\n")
cat("Medium vs Low risk:", exp(summary_fg_muerte$coef[1,1]), "\n")
cat("High vs Low risk:", exp(summary_fg_muerte$coef[2,1]), "\n")

# ------------------------------------------------------------------------------
# 12.2 Discusión de limitaciones del modelo
# ------------------------------------------------------------------------------

cat("\n12.2 Limitaciones del modelo Fine-Gray:\n")
cat("- Asume proporcionalidad de los riesgos de subdistribución\n")
cat("- No captura posibles efectos temporales (no proporcionales)\n")
cat("- La interpretación de los coeficientes es menos intuitiva que en Cox\n")
cat("- Cada modelo se ajusta para un solo evento, ignorando correlación entre eventos\n")
cat("- Sensible al supuesto de independencia condicional entre eventos\n")

# ------------------------------------------------------------------------------
# 12.3 Gráfico de funciones de incidencia acumulada (CIF)
# ------------------------------------------------------------------------------

# Predecimos las CIF para cada grupo de score
cif_recaida_low <- predict(fg_recaida, matrix(c(0,0), nrow = 1))
cif_recaida_medium <- predict(fg_recaida, matrix(c(1,0), nrow = 1))
cif_recaida_high <- predict(fg_recaida, matrix(c(0,1), nrow = 1))

cif_muerte_low <- predict(fg_muerte, matrix(c(0,0), nrow = 1))
cif_muerte_medium <- predict(fg_muerte, matrix(c(1,0), nrow = 1))
cif_muerte_high <- predict(fg_muerte, matrix(c(0,1), nrow = 1))

par(mfrow=c(1,2))

plot(cif_recaida_low, ylim = c(0, 1), col="green", xlab = "Tiempo", ylab = "Probabilidad de sufrir recaída", main = "CIFs según Score")
lines(cif_recaida_medium, lty=2, col="blue")
lines(cif_recaida_high, lty=3, col="red")
legend("topright", legend = c("Low", "Medium", "High"),
       lty = 1:3, col=c("green", "blue", "red"))


plot(cif_muerte_low, ylim = c(0, 1), col="green", xlab = "Tiempo", ylab = "Probabilidad de muerte en remisión", main = "CIFs según Score")
lines(cif_muerte_medium, lty=2, col="blue")
lines(cif_muerte_high, lty=3, col="red")
legend("topright", legend = c("Low", "Medium", "High"),
       lty = 1:3, col=c("green", "blue", "red"))


# ------------------------------------------------------------------------------
# 12.4 Gráfico de bondad de ajuste
# ------------------------------------------------------------------------------

# Usamos el enfoque de residuos de Schoenfeld ponderados
# Para el evento de recaída (similar se puede hacer para muerte)

# Convertimos a formato Surv
datos_surv <- Surv(datos$time, datos$status == 1)

# Ajustamos modelo de Cox con interacción tiempo-score para evaluar no proporcionalidad
cox_test <- coxph(datos_surv ~ score + tt(score), 
                  data = datos,
                  tt = function(x, t, ...) x * log(t))

# Graficamos los residuos de Schoenfeld
ggcoxzph(cox.zph(coxph(datos_surv ~ score, data = datos)))

