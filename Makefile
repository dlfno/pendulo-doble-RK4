# --- Variables de Configuración ---
FC = gfortran
# Flags: -O3 para optimización máxima, -Wall para advertencias
FFLAGS = -O3 -Wall 
TARGET = pendulo
SRC = src/PD_Alonso_CF.f90
DAT = doble_pendulo.dat
VIDEO = pendulo_doble.mp4

# --- Reglas Principales ---

# Regla por defecto: Compilar
all: $(TARGET)

# Compilación del ejecutable
$(TARGET): $(SRC)
	$(FC) $(FFLAGS) -o $(TARGET) $(SRC)

# Ejecutar la simulación numérca
run: $(TARGET)
	./$(TARGET)
	@echo "------------------------------------------------"
	@echo "Simulación completada."
	@echo "Datos guardados en: $(DAT)"
	@echo "------------------------------------------------"

# Visualizar la trayectoria (Caos)
plot: $(DAT)
	gnuplot -p plots/plot_trayectoria.gp

# Generar video MP4 (requiere Python + FFmpeg)
video: $(DAT)
	@echo "Generando animación..."
	python3 src/animate.py
	@echo "Video generado: $(VIDEO)"

# --- VALIDACIÓN FÍSICA ---
# Verifica la conservación de la energía
check: $(DAT)
	@echo "Abriendo gráfico de energía..."
	@echo "Verifica que la línea sea horizontal (constante)."
	gnuplot -p plots/plot_energia.gp

# Limpieza del proyecto
clean:
	rm -f $(TARGET) *.dat *.o *.mod *.mp4

.PHONY: all run plot video check clean
