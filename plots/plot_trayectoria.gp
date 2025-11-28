# Configuración física (deben coincidir con tu Fortran)
L1 = 1.0
L2 = 1.0

# Funciones de transformación de coordenadas
# Masa 1
x1(t1) = L1 * sin(t1)
y1(t1) = -L1 * cos(t1)

# Masa 2 (depende de la posición de la 1)
x2(t1, t2) = x1(t1) + L2 * sin(t2)
y2(t1, t2) = y1(t1) - L2 * cos(t2)

# Configuración del gráfico
set title "Trayectoria de la Masa 2 (Espacio Real)"
set xlabel "Posición X (m)"
set ylabel "Posición Y (m)"
set grid
set size ratio -1  # Para que los ejes tengan la misma escala (no se deforme el círculo)

# Graficamos usando las columnas 2 (t1) y 3 (t2) como argumentos para las funciones x2 y y2
plot "doble_pendulo.dat" using (x2($2, $3)):(y2($2, $3)) with lines \
     title "Trayectoria M2" lc rgb "blue" lw 1

pause -1 "Presiona Enter para salir"
