!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROGRAMA PENDULO_DOBLE
! Solver RK4 de Doble Precisión
!
PROGRAM PENDULO_DOBLE
 IMPLICIT NONE
 
 ! --- Definición de Doble Precisión ---
 INTEGER, PARAMETER :: DP = KIND(1.0D0) 

 INTEGER, PARAMETER :: N = 20001 ! Aumentamos pasos para aprovechar la precisión
 REAL(DP), PARAMETER :: G=9.81_DP, L1=1.0_DP, L2=1.0_DP, M1=1.0_DP, M2=1.0_DP
 REAL(DP) :: H, T, T_MAX
 
 ! Variables temporales para RK4
 REAL(DP) :: DK1_1, DK1_2, DK1_3, DK1_4
 REAL(DP) :: DK2_1, DK2_2, DK2_3, DK2_4
 REAL(DP) :: DK3_1, DK3_2, DK3_3, DK3_4
 REAL(DP) :: DK4_1, DK4_2, DK4_3, DK4_4
 
 ! Estado temporal y Energía
 REAL(DP) :: Y1_TMP, Y2_TMP, Y3_TMP, Y4_TMP, E_NOW
 
 ! Vector de estado: 1=th1, 2=w1, 3=th2, 4=w2
 REAL(DP), DIMENSION (4,N) :: Y 
 INTEGER :: I

 ! --- Interfaces ---
 INTERFACE
    FUNCTION F_THETA(OMEGA)
       IMPORT DP
       REAL(DP) :: F_THETA, OMEGA
    END FUNCTION F_THETA

    FUNCTION F_OMEGA1(TH1, W1, TH2, W2, G, L1, L2, M1, M2)
       IMPORT DP
       REAL(DP) :: F_OMEGA1, TH1, W1, TH2, W2, G, L1, L2, M1, M2
    END FUNCTION F_OMEGA1

    FUNCTION F_OMEGA2(TH1, W1, TH2, W2, G, L1, L2, M1, M2)
       IMPORT DP
       REAL(DP) :: F_OMEGA2, TH1, W1, TH2, W2, G, L1, L2, M1, M2
    END FUNCTION F_OMEGA2

    FUNCTION CALC_ENERGY(TH1, W1, TH2, W2, G, L1, L2, M1, M2)
       IMPORT DP
       REAL(DP) :: CALC_ENERGY, TH1, W1, TH2, W2, G, L1, L2, M1, M2
    END FUNCTION CALC_ENERGY
 END INTERFACE

 ! --- Inicialización ---
 T_MAX = 20.0_DP
 H = T_MAX / REAL(N-1, DP)
 
 ! Condiciones Iniciales (Péndulo casi horizontal)
 Y(1,1) = 1.57_DP ! 90 grados
 Y(2,1) = 0.0_DP
 Y(3,1) = 1.57_DP
 Y(4,1) = 0.0_DP

 ! --- Bucle Runge-Kutta 4 ---
 DO I = 1, N-1
    T = REAL(I-1, DP)*H
    
    Y1_TMP = Y(1,I)
    Y2_TMP = Y(2,I)
    Y3_TMP = Y(3,I)
    Y4_TMP = Y(4,I)

    ! -- K1 --
    DK1_1 = H * F_THETA(Y2_TMP)
    DK1_2 = H * F_OMEGA1(Y1_TMP, Y2_TMP, Y3_TMP, Y4_TMP, G, L1, L2, M1, M2)
    DK1_3 = H * F_THETA(Y4_TMP)
    DK1_4 = H * F_OMEGA2(Y1_TMP, Y2_TMP, Y3_TMP, Y4_TMP, G, L1, L2, M1, M2)

    ! -- K2 --
    DK2_1 = H * F_THETA(Y2_TMP + DK1_2/2.0_DP)
    DK2_2 = H * F_OMEGA1(Y1_TMP + DK1_1/2.0_DP, Y2_TMP + DK1_2/2.0_DP, &
                         Y3_TMP + DK1_3/2.0_DP, Y4_TMP + DK1_4/2.0_DP, G, L1, L2, M1, M2)
    DK2_3 = H * F_THETA(Y4_TMP + DK1_4/2.0_DP)
    DK2_4 = H * F_OMEGA2(Y1_TMP + DK1_1/2.0_DP, Y2_TMP + DK1_2/2.0_DP, &
                         Y3_TMP + DK1_3/2.0_DP, Y4_TMP + DK1_4/2.0_DP, G, L1, L2, M1, M2)

    ! -- K3 --
    DK3_1 = H * F_THETA(Y2_TMP + DK2_2/2.0_DP)
    DK3_2 = H * F_OMEGA1(Y1_TMP + DK2_1/2.0_DP, Y2_TMP + DK2_2/2.0_DP, &
                         Y3_TMP + DK2_3/2.0_DP, Y4_TMP + DK2_4/2.0_DP, G, L1, L2, M1, M2)
    DK3_3 = H * F_THETA(Y4_TMP + DK2_4/2.0_DP)
    DK3_4 = H * F_OMEGA2(Y1_TMP + DK2_1/2.0_DP, Y2_TMP + DK2_2/2.0_DP, &
                         Y3_TMP + DK2_3/2.0_DP, Y4_TMP + DK2_4/2.0_DP, G, L1, L2, M1, M2)

    ! -- K4 --
    DK4_1 = H * F_THETA(Y2_TMP + DK3_2)
    DK4_2 = H * F_OMEGA1(Y1_TMP + DK3_1, Y2_TMP + DK3_2, &
                         Y3_TMP + DK3_3, Y4_TMP + DK3_4, G, L1, L2, M1, M2)
    DK4_3 = H * F_THETA(Y4_TMP + DK3_4)
    DK4_4 = H * F_OMEGA2(Y1_TMP + DK3_1, Y2_TMP + DK3_2, &
                         Y3_TMP + DK3_3, Y4_TMP + DK3_4, G, L1, L2, M1, M2)

    ! -- Actualización --
    Y(1,I+1) = Y(1,I) + (DK1_1 + 2.0_DP*DK2_1 + 2.0_DP*DK3_1 + DK4_1) / 6.0_DP
    Y(2,I+1) = Y(2,I) + (DK1_2 + 2.0_DP*DK2_2 + 2.0_DP*DK3_2 + DK4_2) / 6.0_DP
    Y(3,I+1) = Y(3,I) + (DK1_3 + 2.0_DP*DK2_3 + 2.0_DP*DK3_3 + DK4_3) / 6.0_DP
    Y(4,I+1) = Y(4,I) + (DK1_4 + 2.0_DP*DK2_4 + 2.0_DP*DK3_4 + DK4_4) / 6.0_DP

 END DO

 ! --- Escritura (Formato ES para notación científica precisa) ---
 OPEN(UNIT=10, FILE='doble_pendulo.dat')
 DO I = 1, N
    E_NOW = CALC_ENERGY(Y(1,I), Y(2,I), Y(3,I), Y(4,I), G, L1, L2, M1, M2)
    ! Usamos 18 decimales de precisión en la salida
    WRITE (10,"(4ES25.16)") REAL(I-1, DP)*H, Y(1,I), Y(3,I), E_NOW
 END DO
 CLOSE(10)
 PRINT *, "Simulacion en Doble Precision terminada."

END PROGRAM PENDULO_DOBLE

! -----------------------------------------------------------------
! FUNCIONES (Todas actualizadas a REAL(DP))
! -----------------------------------------------------------------

FUNCTION F_THETA(OMEGA) RESULT(DTH)
 IMPLICIT NONE
 INTEGER, PARAMETER :: DP = KIND(1.0D0)
 REAL(DP) :: OMEGA, DTH
 DTH = OMEGA
END FUNCTION F_THETA

FUNCTION F_OMEGA1(TH1, W1, TH2, W2, G, L1, L2, M1, M2) RESULT(DW1)
 IMPLICIT NONE
 INTEGER, PARAMETER :: DP = KIND(1.0D0)
 REAL(DP) :: TH1, W1, TH2, W2, G, L1, L2, M1, M2, DW1
 REAL(DP) :: DELTA, NUM, DEN
 
 DELTA = TH1 - TH2
 NUM = -G*(2.0_DP*M1+M2)*SIN(TH1) - M2*G*SIN(TH1-2.0_DP*TH2) &
       - 2.0_DP*SIN(DELTA)*M2*(W2*W2*L2 + W1*W1*L1*COS(DELTA))
 DEN = L1 * (2.0_DP*M1 + M2 - M2*COS(2.0_DP*DELTA))
 
 DW1 = NUM / DEN
END FUNCTION F_OMEGA1

FUNCTION F_OMEGA2(TH1, W1, TH2, W2, G, L1, L2, M1, M2) RESULT(DW2)
 IMPLICIT NONE
 INTEGER, PARAMETER :: DP = KIND(1.0D0)
 REAL(DP) :: TH1, W1, TH2, W2, G, L1, L2, M1, M2, DW2
 REAL(DP) :: DELTA, NUM, DEN

 DELTA = TH1 - TH2
 NUM = 2.0_DP*SIN(DELTA) * (W1*W1*L1*(M1+M2) + G*(M1+M2)*COS(TH1) &
       + W2*W2*L2*M2*COS(DELTA))
 DEN = L2 * (2.0_DP*M1 + M2 - M2*COS(2.0_DP*DELTA))

 DW2 = NUM / DEN
END FUNCTION F_OMEGA2

FUNCTION CALC_ENERGY(TH1, W1, TH2, W2, G, L1, L2, M1, M2) RESULT(E_TOTAL)
 IMPLICIT NONE
 INTEGER, PARAMETER :: DP = KIND(1.0D0)
 REAL(DP) :: TH1, W1, TH2, W2, G, L1, L2, M1, M2, E_TOTAL
 REAL(DP) :: T_KIN, V_POT

 T_KIN = 0.5_DP*(M1+M2)*L1*L1*W1*W1 + 0.5_DP*M2*L2*L2*W2*W2 + &
         M2*L1*L2*W1*W2*COS(TH1-TH2)
 
 V_POT = -(M1+M2)*G*L1*COS(TH1) - M2*G*L2*COS(TH2)

 E_TOTAL = T_KIN + V_POT
END FUNCTION CALC_ENERGY
