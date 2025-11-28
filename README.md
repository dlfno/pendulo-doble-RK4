
# Péndulo Doble - Fortran 90 Solucionador RK4

Simulación numérica de un Péndulo Doble caótico (Sistema Dinámico No Lineal) utilizando el método de Runge-Kutta de 4to orden (RK4).

Este proyecto no solo resuelve las ecuaciones de movimiento, sino que incluye herramientas de validación física (conservación de energía), visualización de trayectorias caóticas y generación de animaciones.

## Características

* **Solver:** Implementación de RK4 en Fortran 90.
* **Validación Física:** Cálculo en tiempo real de la Energía Total (Hamiltoniano) para verificar la precisión numérica.
* **Visualización Estática:** Scripts de Gnuplot para ver la evolución de ángulos y la trayectoria en el espacio de fase.
* **Animación:** Script de Python (`matplotlib` + `ffmpeg`) para renderizar video MP4 de la simulación.
* **Automatización:** Gestión completa del flujo de trabajo mediante `Makefile`.

## Requisitos

Para utilizar todas las funciones del repositorio necesitas:

* **Fortran Compiler:** `gfortran` (GCC).
* **Visualización:** `gnuplot`.
* **Automatización:** `make`.
* **Animación (Opcional):**
    * Python 3 (`numpy`, `matplotlib`).
    * FFmpeg (para guardar el video).

## Estructura del Proyecto

```text
.
├── Makefile             # Automatización de compilación y tareas
├── README.md            # Documentación
├── src/
│   ├── PD_Alonso_CF.f90 # Código fuente (Solver RK4 + Física)
│   └── animate.py       # Generador de video Python
└── plots/
    ├── plot_trayectoria.gp  # Grafica el "dibujo" del caos (XY)
    └── plot_energia.gp      # Grafica la estabilidad de la energía
````

## Teoría Matemática

El sistema se modela utilizando la mecánica Lagrangiana. El estado del sistema está definido por los ángulos $\theta_1$ y $\theta_2$.

### Lagrangiano

La energía cinética ($T$) y potencial ($V$) del sistema están dadas por:

$$
T = \frac{1}{2}(m_1 + m_2)l_1^2 \dot{\theta}_1^2 + \frac{1}{2}m_2 l_2^2 \dot{\theta}_2^2 + m_2 l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos(\theta_1 - \theta_2)
$$

$$
V = -(m_1 + m_2)g l_1 \cos \theta_1 - m_2 g l_2 \cos \theta_2
$$

### Ecuaciones de Movimiento

Se resuelven numéricamente las ecuaciones de Euler-Lagrange para las aceleraciones angulares $\ddot{\theta}_1$ y $\ddot{\theta}_2$ (donde $\delta = \theta_1 - \theta_2$):

$$
\ddot{\theta}_1 = \frac{-g(2m_1+m_2)\sin\theta_1 - m_2 g \sin(\theta_1-2\theta_2) - 2\sin\delta \cdot m_2(\dot{\theta}_2^2 l_2 + \dot{\theta}_1^2 l_1 \cos\delta)}{l_1(2m_1+m_2 - m_2\cos(2\delta))}
$$

$$
\ddot{\theta}_2 = \frac{2\sin\delta \cdot (\dot{\theta}_1^2 l_1(m_1+m_2) + g(m_1+m_2)\cos\theta_1 + \dot{\theta}_2^2 l_2 m_2 \cos\delta)}{l_2(2m_1+m_2 - m_2\cos(2\delta))}
$$

## Uso e Instrucciones

El proyecto se gestiona completamente a través del `Makefile`.

### 1\. Compilar y Ejecutar

Compila el código Fortran y ejecuta la simulación numérca. Esto genera el archivo de datos `doble_pendulo.dat`.

```bash
make run
```

### 2\. Validación de Energía (¡Importante\!)

Verifica la calidad de la simulación. Al ser un sistema conservativo (sin fricción), la energía total debe permanecer constante.

```bash
make check
```

  * **Resultado esperado:** Una línea horizontal.
  * **Si la línea sube/baja:** El paso de tiempo `H` es demasiado grande; la simulación está introduciendo errores numéricos.

### 3\. Visualizar Trayectoria

Muestra el gráfico de la trayectoria de la segunda masa en el plano XY.

```bash
make plot
```

### 4\. Generar Video (Animación)

Crea un archivo `pendulo_doble.mp4` visualizando el movimiento.

```bash
make video
```

### 5\. Limpieza

Elimina ejecutables, archivos de datos y videos generados.

```bash
make clean
```

## Formato de Salida (`.dat`)

El archivo `doble_pendulo.dat` contiene 4 columnas:

1.  Tiempo ($t$)
2.  Ángulo 1 ($\theta_1$)
3.  Ángulo 2 ($\theta_2$)
4.  Energía Total ($E$)

-----

*Desarrollado como adaptación de algoritmos de disparo a problemas de valor inicial (IVP) en Física Computacional.*

**Autor:** Alonso Delfino Cervantes Flores
**Licencia:** MIT
