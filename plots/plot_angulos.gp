# Configuración general
set title "Péndulo Doble: Evolución de Ángulos"
set xlabel "Tiempo (s)"
set ylabel "Ángulo (rad)"
set grid

# Estilo de líneas
set style data lines

# Graficar Theta1 (columna 2) y Theta2 (columna 3) vs Tiempo (columna 1)
plot "doble_pendulo.dat" using 1:2 title "Theta 1 (Interno)" lw 2, \
     "doble_pendulo.dat" using 1:3 title "Theta 2 (Externo)" lw 1 lc rgb "red"

pause -1 "Presiona Enter para salir"
