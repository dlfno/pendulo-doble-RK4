!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROGRAMA PENDULO_DOBLE
! Basado en disparo.f90 pero adaptado para un sistema de 4 variables (IVP)
!
PROGRAM PENDULO_DOBLE
 IMPLICIT NONE
 INTEGER, PARAMETER :: N=10001 ! Más pasos para mejor resolución temporal
 REAL, PARAMETER :: G=9.81, L1=1.0, L2=1.0, M1=1.0, M2=1.0
 REAL :: H, T, T_MAX
 ! Variables temporales para RK4 (4 variables x 4 pasos k)
 REAL :: DK1_1, DK1_2, DK1_3, DK1_4
 REAL :: DK2_1, DK2_2, DK2_3, DK2_4
 REAL :: DK3_1, DK3_2, DK3_3, DK3_4
 REAL :: DK4_1, DK4_2, DK4_3, DK4_4
 ! Estado temporal
 REAL :: Y1_TMP, Y2_TMP, Y3_TMP, Y4_TMP
 ! Vector de estado: Indices: 1=theta1, 2=omega1, 3=theta2, 4=omega2
 REAL, DIMENSION (4,N) :: Y 
 INTEGER :: I
 REAL :: E_NOW ! Variable para almacenar la energía actual

 ! --- Funciones de derivadas ---
 INTERFACE
    FUNCTION F_THETA(OMEGA)
       REAL :: F_THETA, OMEGA
    END FUNCTION F_THETA
    FUNCTION F_OMEGA1(TH1, W1, TH2, W2, G, L1, L2, M1, M2)
       REAL :: F_OMEGA1, TH1, W1, TH2, W2, G, L1, L2, M1, M2
    END FUNCTION F_OMEGA1
    FUNCTION F_OMEGA2(TH1, W1, TH2, W2, G, L1, L2, M1, M2)
       REAL :: F_OMEGA2, TH1, W1, TH2, W2, G, L1, L2, M1, M2
    END FUNCTION F_OMEGA2
 ! --- Energía ---
    FUNCTION CALC_ENERGY(TH1, W1, TH2, W2, G, L1, L2, M1, M2)
       REAL :: CALC_ENERGY, TH1, W1, TH2, W2, G, L1, L2, M1, M2
    END FUNCTION CALC_ENERGY
 END INTERFACE

 ! --- Inicialización ---
 T_MAX = 20.0       ! Tiempo total de simulación
 H = T_MAX / (N-1)  ! Paso de tiempo (dt)
 
 ! Condiciones Iniciales (en radianes y rad/s)
 Y(1,1) = 1.57 ! Theta1 inicial (90 grados)
 Y(2,1) = 0.0  ! Omega1 inicial
 Y(3,1) = 1.57 ! Theta2 inicial (90 grados)
 Y(4,1) = 0.0  ! Omega2 inicial

 ! --- Bucle Runge-Kutta de 4to Orden ---
 DO I = 1, N-1
    T = (I-1)*H
    
    ! Valores actuales
    Y1_TMP = Y(1,I)
    Y2_TMP = Y(2,I)
    Y3_TMP = Y(3,I)
    Y4_TMP = Y(4,I)

    ! -- K1 --
    DK1_1 = H * F_THETA(Y2_TMP)
    DK1_2 = H * F_OMEGA1(Y1_TMP, Y2_TMP, Y3_TMP, Y4_TMP, G, L1, L2, M1, M2)
    DK1_3 = H * F_THETA(Y4_TMP)
    DK1_4 = H * F_OMEGA2(Y1_TMP, Y2_TMP, Y3_TMP, Y4_TMP, G, L1, L2, M1, M2)

    ! -- K2 -- (Usando T + H/2 y Y + K1/2)
    DK2_1 = H * F_THETA(Y2_TMP + DK1_2/2.0)
    DK2_2 = H * F_OMEGA1(Y1_TMP + DK1_1/2.0, Y2_TMP + DK1_2/2.0, &
                         Y3_TMP + DK1_3/2.0, Y4_TMP + DK1_4/2.0, G, L1, L2, M1, M2)
    DK2_3 = H * F_THETA(Y4_TMP + DK1_4/2.0)
    DK2_4 = H * F_OMEGA2(Y1_TMP + DK1_1/2.0, Y2_TMP + DK1_2/2.0, &
                         Y3_TMP + DK1_3/2.0, Y4_TMP + DK1_4/2.0, G, L1, L2, M1, M2)

    ! -- K3 -- (Usando T + H/2 y Y + K2/2)
    DK3_1 = H * F_THETA(Y2_TMP + DK2_2/2.0)
    DK3_2 = H * F_OMEGA1(Y1_TMP + DK2_1/2.0, Y2_TMP + DK2_2/2.0, &
                         Y3_TMP + DK2_3/2.0, Y4_TMP + DK2_4/2.0, G, L1, L2, M1, M2)
    DK3_3 = H * F_THETA(Y4_TMP + DK2_4/2.0)
    DK3_4 = H * F_OMEGA2(Y1_TMP + DK2_1/2.0, Y2_TMP + DK2_2/2.0, &
                         Y3_TMP + DK2_3/2.0, Y4_TMP + DK2_4/2.0, G, L1, L2, M1, M2)

    ! -- K4 -- (Usando T + H y Y + K3)
    DK4_1 = H * F_THETA(Y2_TMP + DK3_2)
    DK4_2 = H * F_OMEGA1(Y1_TMP + DK3_1, Y2_TMP + DK3_2, &
                         Y3_TMP + DK3_3, Y4_TMP + DK3_4, G, L1, L2, M1, M2)
    DK4_3 = H * F_THETA(Y4_TMP + DK3_4)
    DK4_4 = H * F_OMEGA2(Y1_TMP + DK3_1, Y2_TMP + DK3_2, &
                         Y3_TMP + DK3_3, Y4_TMP + DK3_4, G, L1, L2, M1, M2)

    ! -- Actualización de Estado --
    Y(1,I+1) = Y(1,I) + (DK1_1 + 2.0*DK2_1 + 2.0*DK3_1 + DK4_1) / 6.0
    Y(2,I+1) = Y(2,I) + (DK1_2 + 2.0*DK2_2 + 2.0*DK3_2 + DK4_2) / 6.0
    Y(3,I+1) = Y(3,I) + (DK1_3 + 2.0*DK2_3 + 2.0*DK3_3 + DK4_3) / 6.0
    Y(4,I+1) = Y(4,I) + (DK1_4 + 2.0*DK2_4 + 2.0*DK3_4 + DK4_4) / 6.0

 END DO

! --- Escritura de resultados ---
 OPEN(UNIT=10, FILE='doble_pendulo.dat')
 DO I = 1, N
    ! Calculamos la energía en este paso
    E_NOW = CALC_ENERGY(Y(1,I), Y(2,I), Y(3,I), Y(4,I), G, L1, L2, M1, M2)
    
    ! Escribe: Tiempo, Theta1, Theta2, Energia
    WRITE (10,"(4F16.8)") (I-1)*H, Y(1,I), Y(3,I), E_NOW
 END DO
 CLOSE(10)
 PRINT *, "Calculo terminado. Datos (incluyendo energia) en 'doble_pendulo.dat'"

END PROGRAM PENDULO_DOBLE

! -----------------------------------------------------------------
! FUNCIONES DEL SISTEMA FISICO
! -----------------------------------------------------------------

! Derivada de Theta es Omega
FUNCTION F_THETA(OMEGA) RESULT(DTH)
 IMPLICIT NONE
 REAL :: OMEGA, DTH
 DTH = OMEGA
END FUNCTION F_THETA

! Ecuación de movimiento para Omega1 (Aceleración angular 1)
FUNCTION F_OMEGA1(TH1, W1, TH2, W2, G, L1, L2, M1, M2) RESULT(DW1)
 IMPLICIT NONE
 REAL :: TH1, W1, TH2, W2, G, L1, L2, M1, M2, DW1
 REAL :: DELTA, NUM, DEN
 
 DELTA = TH1 - TH2
 ! Ecuación simplificada del Lagrangiano
 NUM = -G*(2.0*M1+M2)*SIN(TH1) - M2*G*SIN(TH1-2.0*TH2) &
       - 2.0*SIN(DELTA)*M2*(W2*W2*L2 + W1*W1*L1*COS(DELTA))
 DEN = L1 * (2.0*M1 + M2 - M2*COS(2.0*DELTA))
 
 DW1 = NUM / DEN
END FUNCTION F_OMEGA1

! Ecuación de movimiento para Omega2 (Aceleración angular 2)
FUNCTION F_OMEGA2(TH1, W1, TH2, W2, G, L1, L2, M1, M2) RESULT(DW2)
 IMPLICIT NONE
 REAL :: TH1, W1, TH2, W2, G, L1, L2, M1, M2, DW2
 REAL :: DELTA, NUM, DEN

 DELTA = TH1 - TH2
 NUM = 2.0*SIN(DELTA) * (W1*W1*L1*(M1+M2) + G*(M1+M2)*COS(TH1) &
       + W2*W2*L2*M2*COS(DELTA))
 DEN = L2 * (2.0*M1 + M2 - M2*COS(2.0*DELTA))

 DW2 = NUM / DEN
END FUNCTION F_OMEGA2

! Calcula la Energía Total (Cinética + Potencial)
FUNCTION CALC_ENERGY(TH1, W1, TH2, W2, G, L1, L2, M1, M2) RESULT(E_TOTAL)
 IMPLICIT NONE
 REAL :: TH1, W1, TH2, W2, G, L1, L2, M1, M2, E_TOTAL
 REAL :: T_KIN, V_POT

 ! Energía Cinética (T)
 T_KIN = 0.5*(M1+M2)*L1*L1*W1*W1 + 0.5*M2*L2*L2*W2*W2 + &
         M2*L1*L2*W1*W2*COS(TH1-TH2)
 
 ! Energía Potencial (V) - Asumiendo y=0 en el pivote y eje Y hacia arriba
 V_POT = -(M1+M2)*G*L1*COS(TH1) - M2*G*L2*COS(TH2)

 E_TOTAL = T_KIN + V_POT
END FUNCTION CALC_ENERGY
