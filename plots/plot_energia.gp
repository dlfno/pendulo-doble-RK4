# Configuración
set title "Conservación de la Energía (Hamiltoniano)"
set xlabel "Tiempo (s)"
set ylabel "Energía Total (J)"
set grid

# Formato de números en ejes
set format y "%.4f"

# Estadísticas para centrar la gráfica (opcional)
stats "doble_pendulo.dat" using 4 nooutput
E_mean = STATS_mean

# Graficar la Energía (Columna 4)
# Usamos un rango ajustado para ver las pequeñas variaciones
plot "doble_pendulo.dat" using 1:4 title "Energía Total" lw 2 lc rgb "green"

pause -1 "Presiona Enter para salir"
