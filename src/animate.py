import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- Configuración ---
ARCHIVO_DATOS = 'doble_pendulo.dat'
L1 = 1.0  # Longitud varilla 1 (debe coincidir con Fortran)
L2 = 1.0  # Longitud varilla 2
FPS = 30  # Cuadros por segundo del video
SKIP = 5  # Saltar pasos para acelerar la renderización (reduce tiempo de generación)

print(f"Leyendo datos de {ARCHIVO_DATOS}...")
try:
    data = np.loadtxt(ARCHIVO_DATOS)
except OSError:
    print("Error: No se encuentra el archivo .dat. Ejecuta primero el código Fortran.")
    exit()

# Extraer columnas (Tiempo, Theta1, Theta2)
# data[:, 0] es tiempo
th1 = data[::SKIP, 1] 
th2 = data[::SKIP, 2]

# --- Transformación a Coordenadas Cartesianas ---
# Posición de la masa 1
x1 = L1 * np.sin(th1)
y1 = -L1 * np.cos(th1)

# Posición de la masa 2
x2 = x1 + L2 * np.sin(th2)
y2 = y1 - L2 * np.cos(th2)

# --- Configuración de la Figura ---
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-(L1+L2)-0.5, (L1+L2)+0.5)
ax.set_ylim(-(L1+L2)-0.5, (L1+L2)+0.5)
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)
ax.set_title(f'Simulación Péndulo Doble (RK4)')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')

# Elementos gráficos a animar
line, = ax.plot([], [], 'o-', lw=2, color='black') # Las varillas
trace, = ax.plot([], [], '-', lw=1, color='red', alpha=0.5) # La estela
time_template = 'Tiempo = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

# Historial para la estela (últimos 100 puntos)
history_len = 100
history_x, history_y = [], []

def init():
    line.set_data([], [])
    trace.set_data([], [])
    time_text.set_text('')
    return line, trace, time_text

def update(i):
    # Coordenadas actuales de las masas (origen -> m1 -> m2)
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    # Actualizar varillas
    line.set_data(thisx, thisy)

    # Actualizar estela (trail)
    if i == 0:
        history_x.clear()
        history_y.clear()
    
    history_x.append(x2[i])
    history_y.append(y2[i])
    
    # Mantener solo los últimos N puntos para no saturar
    if len(history_x) > history_len:
        history_x.pop(0)
        history_y.pop(0)
        
    trace.set_data(history_x, history_y)
    
    # Actualizar tiempo (aproximado basado en frames)
    # Nota: Si quieres tiempo exacto, extrae la columna 0 de datos
    current_time = data[i*SKIP, 0]
    time_text.set_text(time_template % current_time)
    
    return line, trace, time_text

# Crear animación
print(f"Generando animación ({len(th1)} cuadros)...")
ani = animation.FuncAnimation(
    fig, update, frames=len(th1),
    init_func=init, interval=1000/FPS, blit=True
)

# Guardar como MP4
output_file = 'pendulo_doble.mp4'
print(f"Guardando {output_file}...")
ani.save(output_file, writer='ffmpeg', fps=FPS, dpi=150)
print("¡Listo!")
