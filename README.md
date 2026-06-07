# Numeric-Arm-Figure-Tracer

Simulación de un brazo robótico de 3 grados de libertad (3R) capaz de seguir trayectorias extraídas de imágenes mediante procesamiento digital de señales y cinemática inversa.

## Estructura del Proyecto
- `trazos_por_parte.m`: Procesa una imagen, detecta bordes (Canny) y genera las trayectorias en formato CSV.
- `main_simulacion.m`: Script principal que carga la trayectoria, resuelve la cinemática inversa y ejecuta la animación.
- `forward_kinematics.m`: Cálculo de la posición del efector final a partir de los ángulos de las juntas.
- `inverse_kinematics_solver.m`: Implementación del método de la Pseudoinversa de Moore-Penrose para resolver la cinemática inversa.
- `jacobian_matrix.m`: Cálculo de la matriz Jacobiana del sistema.

## Teoría y Métodos
### Cinemática Inversa
El brazo es redundante (3 juntas para 2 dimensiones). Se utiliza la matriz Jacobiana $J(\theta)$ para relacionar las velocidades de las juntas con las velocidades del efector final:
$$\dot{x} = J(\theta) \dot{\theta}$$
Para hallar los ángulos, se aplica el método iterativo basado en la pseudoinversa:
$$\Delta\theta = J^\dagger \Delta x$$

### Procesamiento de Imagen
Se utiliza el detector de bordes de Canny para extraer contornos precisos de figuras. Estos contornos se filtran por longitud y se remuestrean para asegurar una velocidad de dibujo constante en la simulación.

## Requisitos
- MATLAB (Probado en R2025b)
- Image Processing Toolbox

## Uso
1. Ejecuta `trazos_por_parte.m` para generar el archivo `trayectoria_all.csv` a partir de una imagen.
2. Ejecuta `main_simulacion.m` para ver la animación del brazo dibujando la figura.
