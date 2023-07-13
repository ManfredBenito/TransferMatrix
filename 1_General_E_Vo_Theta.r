####################################
######## NOTAS IMPORTANTES #########
####################################

# THETA ES EL ÁNGULO DE INCIDENCIA, NO EL DE TRANSMISIÓN
# Es decir, una función que pida el ángulo theta,
# requiere únicamente el ángulo de incidencia para realizar los cálculos


####################################
############ CONSTANTES##############
####################################

########### CONSTANTES #############

vF <- 10**15 # Velocidad de Fermi [nm/s]
hbar <- 6.58 * 10**(-16) # h_barra [eV*s]

###### PARÁMETROS INICIALES ########

# theta <- 0
N <- 1 # Barreras
xo <- 0 # Origen en el eje x
a <- 0.142 # Separación entre barreras [nm]
L <- 0.178 # Ancho de las barreras [nm]
#Vo <- 100 # Potencial (Altura) de las barreras [eV]
#E <- 300 # Energía de la partícula [eV]


### CALCULO DE OTRAS CONSTANTES ###

#### REGIONES #####
R <- N - 1 # Calcula las regiones
# Regiones  2 + 2R

############################################################
######################### FUNCIONES ########################
############################################################

quad <- function(a, b, c) {
  a <- as.complex(a)
  resultado <- # c(
    (-b + sqrt(b^2 - 4 * a * c)) / (2 * a) # ,
  # (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
  # if(all(Im(resultado) == 0)) resultado <- Re(resultado)
  # if(resultado[1] == resultado[2]) return(resultado[1])
  resultado
}

mysqrt <- function(x){
  if (x > 0){
    return( sqrt(x))
  } else {
     return( 1i* sqrt( abs(x)))
  }
}

############################################################
######### SON CONSTANTES DURANTE TODO EL PROGRAMA ##########
############################################################

################## CÁLCULO DE SIGNO ########################

####### CÁLCULO DE POSICIONES GUARDADAS EN UN VECTOR #######

Posiciones <- function(R, a, L, xo) {
  # Esta función regresa las posiciones de donde el potencial cambia
  X <- c(xo, xo + a, xo + a + L)
  if (R < 1) {
    return(X)
  } else {
    for (n in 1:R) {
      # print(n)
      ba <- X[length(X)] + a
      X <- append(X, ba)
      ba <- X[length(X)] + L
      X <- append(X, ba)
    }
  }
  return(X)
}

##########  GENERA LAS BARRERAS DE POTENCIAL ##############
Funcion_barreras <- function(R, a, L, xo, V1, Vb) {
  Y <- c(Vb, V1, Vb)
  if (R < 1) {
    return(Y)
  } else {
    for (n in 1:R) {
      # print(n)
      ba <- V1
      Y <- append(Y, ba)
      ba <- Vb
      Y <- append(Y, ba)
    }
  }
  return(Y)
}

#################################################
############ Generación de posiciones ###########
#################################################

X <- Posiciones(R, a, L, xo)
#print("Posiciones: ")
#print(X)
#print("Longitud del vector de posiciones: ")
#print(length(X))
## print(length(X))
Y <- Funcion_barreras(R, a, L, xo, Vo, 0)
#print("Barreras: ")
#print(Y)
# print(length(Y))

png("Barreras_generadas.png")
plot(X, Y, type = "s", ylab = "E", main = "Gráfica que ilustra la disposición de las barreras")
dev.off()

#############################################################
############## FUNCIONES NO CONSTANTES ######################
#############################################################



####### FUNCIÓN QUE EVALÚA EN LA MATRIZ DE EXPONENCIALES #######
Matriz_evaluada <- function(x, theta = 0, Kx, s = 1) {
  # Crear la matriz A a partir de los parámetros
  A <- matrix(
    c(
      exp(1i * Kx * x),
      s * exp(1i * theta) * exp(1i * Kx * x),
      exp(-1i * Kx * x),
      -s * exp(-1i * theta) * exp(-1i * Kx * x)
    ),
    nrow = 2, ncol = 2
  )
  return(A)
}

####### CALCULO DEL ÁNGULO DE TRANSMISIÓN ####### ¿PUEDE SER CONSTANTE?
Angulo_transmision <- function(theta, Px, Py) {
  angulo_t <- atan(Py/Px)
  #angulo_t <- asin((Ke1 / Ke2) * sin(theta))
  return(angulo_t)
}

####### CÁLCULO DEL COEFICIENTE Kx #######
CoeficienteKx <- function(theta, E, Vo) {
  Ke1 <- (E / (hbar * vF))
  Ke2 <- ((E - Vo) / (hbar * vF))
  #K_x <- Ke1 * cos(theta)
  K_x <- mysqrt(Ke1**2 - (Ke1 * sin(theta))**2)
  return(K_x)
}

####### CÁLCULO DEL COEFICIENTE Qx #######
CoeficienteQx <- function(theta, E, Vo) {
  Ke1 <- (E / (hbar * vF))
  Ke2 <- ((E - Vo) / (hbar * vF))
  Q_x <- mysqrt(Ke2**2 - (Ke1 * sin(theta))**2)
  return(Q_x)
}
####################################################################
###################### Calculadora #################################
####################################################################

Coeficiente_Transmision <- function(Posiciones, theta, E, Vo){
  s_out <- sign(E)
  s_in <- sign( E - Vo)
  
  Kx <- CoeficienteKx(theta, E, Vo)
  Qx <- CoeficienteQx(theta, E, Vo)
  Ky <- (E / (hbar * vF) ) * sin(theta)

  #print(Kx)
  #print(Qx)

  phi <- Angulo_transmision(theta, Qx, Ky)

  A <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  #print(A)
  for (m in 1:(N * 2)) {
    #print("La posición es: ")
    #print(Posiciones[m+1])
    M_out <- Matriz_evaluada(Posiciones[m+1], theta, Kx, s_out)
    M_in <- Matriz_evaluada(Posiciones[m+1], phi, Qx, s_in)

    aux <- m %% 2
    if (aux==1){
      A <- A %*% solve(M_out) %*% M_in
    } else {
      A <- A %*% solve(M_in) %*% M_out
    }
  }

  A_III <- 1 / A[1]

  T <- Mod(A_III) ** 2

  if (T < 10**(-10)){
    T <- 0
  }
  return(T)
}

angulos <- seq(-1 * pi / 2, pi / 2, length = 500)
Y <- c()
for (theta in angulos) {
  T <- Coeficiente_Transmision(X,theta, 75, 100)
  Y <- append(Y, T, length(Y))
  }

png("Coeficiente.png")
plot(angulos, Y, type = "l", main = "Coeficiente de transmisión", xlab = "theta", ylab = "T",ylim=c(-0.2,1))
dev.off()

#print(angulos)

data_frame <- data.frame(
  "Posicińn" = angulos,
  "Coeficiente de Transmisión" = Y
)

write.csv(data_frame, "./TvsX.csv", row.names=TRUE)

Energias <- seq(10, 800, length = 2000)
Y2 <- c()

for (E in Energias){
  T <- Coeficiente_Transmision(X, pi/4, E, 200)
  Y2 <- append(Y2, T, length(Y2))
}

png("CoeficienteVSEnergia.png")
plot(Energias, type='l', Y2, main = "Coeficiente de transmisión", xlab = "E", ylab = "T",ylim=c(-0.2,1))
dev.off()


Energias_T <- data.frame(
  "Energías" = Energias,
  "Coeficiente de Transmisión" = Y2
)

write.csv(Energias_T, "./TvsE.csv", row.names=TRUE)