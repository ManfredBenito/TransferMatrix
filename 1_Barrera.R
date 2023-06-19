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
Vo <- 300 # Potencial (Altura) de las barreras [eV]
E <- 100 # Energía de la partícula [eV]
s <- 1 # ¿Vanda de Conducción o de Valencia?

### CALCULO DE OTRAS CONSTANTES ###

#### REGIONES #####
R <- N - 1 # Calcula las regiones
# Regiones  2 + 2R

### CÁLCULO DE CONSTANTES PARA KX Y QX ###
Ke1 <- (E / (hbar * vF))
Ke2 <- ((E - Vo) / (hbar * vF))

# print(Ke1)
# print(Ke2)


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

############################################################
######### SON CONSTANTES DURANTE TODO EL PROGRAMA ##########
############################################################

################## CÁLCULO DE SIGNO ########################
s1 <- sign(E - Vo)
s2 <- sign(E - 0)

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
print("Posiciones: ")
print(X)
print(length(X))
# print(length(X))
Y <- Funcion_barreras(R, a, L, xo, Vo, 0)
print("Barreras: ")
print(Y)
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
      exp(1i * Kx * x), exp(-1i * Kx * x),
      s * exp(1i * theta) * exp(1i * Kx * x),
      -s * exp(-1i * theta) * exp(-1i * Kx * x)
    ),
    nrow = 2, ncol = 2
  )
  # Devolver la matriz resultante
  # Falta agregar el parametro s que indica en qué bandas estoy trabajando
  return(A)
}

####### CALCULO DEL ÁNGULO DE TRANSMISIÓN ####### ¿PUEDE SER CONSTANTE?
Angulo_transmision <- function(theta) {
  angulo_t <- asin((Ke1 / Ke2) * sin(theta))
  return(angulo_t)
}

####### CÁLCULO DEL COEFICIENTE Kx #######
CoeficienteKx <- function(theta) {
  K_x <- Ke1 * cos(theta)
  return(K_x)
}

####### CÁLCULO DEL COEFICIENTE Qx #######
CoeficienteQx <- function(theta) {
  Q_x <- sqrt(Ke2**2 - (Ke1 * sin(theta))**2)
  return(Q_x)
}

######## ASIGNACIÓN DE LOS VALORES PARA K_x o Q_x #########
Coef_Kx <- function(b, theta, posX = X, posY = Y) { # b, es el 2° índice y theta el ángulo
  aux <- b %% 2
  if (aux == 1) {
    px <- CoeficienteKx(theta)
    # print('Kx')
  } else {
    px <- CoeficienteQx(theta)
    # print('Qx')
  }
  return(px)
}

######### FUNCIÓN DE MATRICES PARA LAS DISTINTAS REGIONES ##########
# Nota: Esta función calcula el signo s = sgn(E-Vo)
# Nota: la función, calcula el ángulo de transmisión dentro de las barreras
# Nota: también calcula y asigna el coeficiente a usar Kx o Qx

Matriz_Region <- function(a, b, theta) {
  x <- X[a]
  # D <- c('a=',a,'b=', b,'x=', x, 'theta=',theta)
  # print(D)

  theta_t <- Angulo_transmision(theta)
  aux <- b %% 2
  if (aux == 1) {
    s <- sign(E)
    angulo <- theta
    # print('Kx')
  } else {
    s <- sign(E - Vo)
    angulo <- theta_t ####### Hacer una función pensando en Generalizar ####
    # print('Qx')
  }

  print("El signo es")
  print(s)

  Px <- Coef_Kx(b, theta)

  # print(Kx)

  if (a != b) { # Revisar el caso de N > 1 M_{2,2} ya no funciona
    # s <- s1
    A <- Matriz_evaluada(x, angulo, Px, s)
    A <- solve(A)
    # print('Se usó Theta:') MAL
    # print(theta_t) MAL
    # print('Sí se invirtió')
  } else {
    # s <- s2
    A <- Matriz_evaluada(x, angulo, Px, s)
    # print('Se usó Theta_t:')
    # print(theta)
  }

  print(A)
  return(A)
}

########## MATRIX IDENTIDAD ##########
# A <- matrix(c(1,0,0,1) ,nrow = 2, ncol = 2)

Funcion_Calculadora <- function(theta, X) {
  A <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  print(A)
  for (i in length(X):2) { #### ERROR ##### En la numeración
    print("Pos evaluada: ")
    print(X[i])
    print("Indice: ")
    print(a)
    M_n_1 <- Matriz_Region(i, i - 1, theta)
    M_n_2 <- Matriz_Region(i, i, theta)
    # write.table(M_n_1,file='Matrices.txt',row.names = F, col.names = F)
    # unlink('Matrices.txt',)
    # write.table(M_n_2,file='Matrices.txt',row.names = F, col.names = F)
    # unlink('Matrices.txt')
    A <- M_n_1 %*% M_n_2 %*% A
    print("La matriz resultante es")
    print(A)
  }
  print("Matriz final")
  print(A)
  print("La solución es  ")
  #  a1 <- A[1]*A[3]
  #  b1 <- -1*A[3]
  #  t <- quad(a1,b1,0)

  t <- 1 / A[1]

  print(t)
  T <- sqrt(t * t)
  # T <- Mod(T)
  T <- Re(T)
  print(T)
  return(T)
}


###################################################

angulos <- seq(-1 * pi / 2, pi / 2, length = 100)
print(typeof(angulos))
# angulos <- c(-1 * pi / 2.5)

Y <- c()

# text('Log.txt')
for (phi in angulos) {
  PMax <- 0
  tparametro <- Funcion_Calculadora(phi, X)
  Y <- append(Y, tparametro, length(Y))
  if (tparametro > PMax) {
    print("PHI MAXIMO")
    print(phi)
  }
  PMax <- phi
}
# dev.off()

png("Coeficiente.png")
plot(angulos, Y, type = "l", main = "Coeficiente de transmisión", xlab = "theta", ylab = "T")
dev.off()

print("Valor máximo")
print(max(Y))
print(PMax)

Observador <- Funcion_Calculadora(0, X)
