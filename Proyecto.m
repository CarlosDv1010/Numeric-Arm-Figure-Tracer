%% SECCIÓN 2: Ejecución de la Simulación

% --- 1. Leer la trayectoria desde el archivo CSV maestro ---
% Asegúrate de que el nombre del archivo coincida con el generado por tu nuevo script.
try
    trayectoria_pixeles = readmatrix('Pichu_trayectoria_all.csv');
catch
    error('No se pudo encontrar el archivo CSV. Ejecuta primero tu nuevo script de extracción.');
end


% --- 2. Escalar y centrar los puntos para el espacio del robot ---
% Esta lógica es clave para adaptar los píxeles a centímetros.

reach_target = 13; % Tamaño objetivo del dibujo en el espacio del robot (en cm)

% Calcular el "bounding box" global, ignorando los NaN
puntos_validos = trayectoria_pixeles(~isnan(trayectoria_pixeles(:,1)), :);
xmin = min(puntos_validos(:,1));
xmax = max(puntos_validos(:,1));
ymin = min(puntos_validos(:,2));
ymax = max(puntos_validos(:,2));

ancho_pix = xmax - xmin;
alto_pix = ymax - ymin;

x_centro_pix = (xmin + xmax) / 2;
y_centro_pix = (ymin + ymax) / 2;

% Calcular el factor de escala
escala = reach_target / max(ancho_pix, alto_pix);

% Aplicar la transformación: centrar, escalar e invertir el eje Y
x_escalado = (trayectoria_pixeles(:,1) - x_centro_pix) * escala;
% Invertimos Y para que la orientación sea correcta (en imágenes Y crece hacia abajo)
y_escalado = -(trayectoria_pixeles(:,2) - y_centro_pix) * escala;


% --- 3. Trasladar al área de trabajo del robot ---
% Movemos el centro del dibujo a una posición cómoda, ej. (12, 0)
x_offset = 0;
y_offset = 0;

x_final = x_escalado + x_offset;
y_final = y_escalado + y_offset;

% Esta es la variable final que usará el resto de tu código
trayectoria_robot_final = [x_final, y_final];


% --- 4. (Opcional pero recomendado) Visualizar la trayectoria escalada ---
figure('Name', 'Trayectoria Final en Espacio del Robot');
plot(trayectoria_robot_final(:,1), trayectoria_robot_final(:,2), 'b.-');
axis equal;
grid on;
title('Trayectoria lista para el Solver de IK');
xlabel('X (cm)');
ylabel('Y (cm)');

fprintf('✅ Trayectoria cargada y escalada correctamente.\n');





% Posición inicial del robot (según el PDF, en grados)
theta_inicial_deg = [10; 0; 0];
theta_actual_rad = deg2rad(theta_inicial_deg); % ¡Siempre trabajar en radianes!

% Almacenaremos todos los ángulos calculados aquí
historial_angulos = {};
% Parámetros del brazo (en cm, como en el PDF)
L1 = 6.0;
L2 = 5.5;
L3 = 3.0;
% Recorrer cada trazo (cada célula en 'trayectorias_robot')
for i = 1:size(trayectoria_robot_final, 1)
    
    punto_actual = trayectoria_robot_final(i, :); % Usar paréntesis ()
    
    % Si la fila actual es un NaN, es un separador de trazos.
    if isnan(punto_actual(1))
        % Guardamos un NaN en el historial de ángulos para la animación
        historial_angulos{end+1} = NaN(3,1);
        continue; % Saltamos al siguiente punto
    end
    
    % Punto objetivo actual [x; y]
    punto_objetivo = [punto_actual(1); punto_actual(2)];
    
    % Resolver la cinemática inversa para este punto   
    % La suposición inicial es la última posición calculada, para suavidad
    theta_actual_rad = inverse_kinematics_solver(punto_objetivo, theta_actual_rad, L1, L2, L3);
    
    % Guardar el resultado
    historial_angulos{end+1} = theta_actual_rad;
end

fprintf('✅ ¡Cálculo de cinemática inversa completado!\n');

% Ahora 'historial_angulos' contiene la secuencia completa de ángulos 
% que tu brazo debe seguir para dibujar la figura.
% ¡Estás listo para la animación y el resto del análisis!



%% ============================================================
%  PARTE 2: SOLVER DE CINEMÁTICA INVERSA Y SIMULACIÓN
% =============================================================

%% SECCIÓN 1: Definiciones del Robot y Funciones Clave

% Función de Cinemática Directa (FK)
% Recibe los 3 ángulos (en radianes) y devuelve la posición [x, y] de la punta.
function pos = forward_kinematics(thetas, L1, L2, L3)
    th1 = thetas(1);
    th2 = thetas(2);
    th3 = thetas(3);
    
    x = L1*cos(th1) + L2*cos(th1+th2) + L3*cos(th1+th2+th3);
    y = L1*sin(th1) + L2*sin(th1+th2) + L3*sin(th1+th2+th3);
    
    pos = [x; y];
end

% Función para la Matriz Jacobiana
% Recibe los 3 ángulos (en radianes) y devuelve la matriz Jacobiana 2x3.
function J = jacobian_matrix(thetas, L1, L2, L3)
    th1 = thetas(1);
    th2 = thetas(2);
    th3 = thetas(3);

    % Derivadas parciales de las ecuaciones de FK
    J11 = -L1*sin(th1) - L2*sin(th1+th2) - L3*sin(th1+th2+th3);
    J12 = -L2*sin(th1+th2) - L3*sin(th1+th2+th3);
    J13 = -L3*sin(th1+th2+th3);

    J21 = L1*cos(th1) + L2*cos(th1+th2) + L3*cos(th1+th2+th3);
    J22 = L2*cos(th1+th2) + L3*cos(th1+th2+th3);
    J23 = L3*cos(th1+th2+th3);
    
    J = [J11, J12, J13; 
         J21, J22, J23];
end

% Función del Solver de Cinemática Inversa (IK)
% Recibe un punto objetivo [x_d; y_d] y los ángulos actuales [th1; th2; th3]
% Devuelve los nuevos ángulos que alcanzan el objetivo.
function thetas_final = inverse_kinematics_solver(target_pos, current_thetas, L1, L2, L3)
    
    thetas = current_thetas; % Empezamos desde la posición actual
    
    % Parámetros del solver
    tolerancia = 0.01;      % Tolerancia de error (en cm)
    max_iteraciones = 100;  % Para evitar bucles infinitos
    
    for i = 1:max_iteraciones
        % 1. Calcular la posición actual del brazo
        current_pos = forward_kinematics(thetas, L1, L2, L3);
        
        % 2. Calcular el vector de error
        error_pos = target_pos - current_pos;
        
        % 3. Verificar si ya llegamos (si el error es muy pequeño)
        if norm(error_pos) < tolerancia
            break; % ¡Éxito! Salimos del bucle
        end
        
        % 4. Calcular el Jacobiano en la posición actual
        J = jacobian_matrix(thetas, L1, L2, L3);
        
        % 5. Calcular la pseudoinversa (¡la clave de Moore-Penrose!)
        J_pinv = pinv(J);
        
        % 6. Calcular el cambio necesario en los ángulos
        delta_thetas = J_pinv * error_pos;
        
        % 7. Actualizar los ángulos
        thetas = thetas + delta_thetas;
    end
    
    thetas_final = thetas;
end


%% ============================================================
%  PARTE 3: ANIMACIÓN Y VISUALIZACIÓN (INICIO DIRECTO)
% =============================================================

% --- PASO 1: Preparar la figura y encontrar el punto de inicio ---

figure('Name', 'Simulación del Brazo Robótico');
hold on; axis equal; grid on;
xlabel('X (cm)'); ylabel('Y (cm)'); title('Figura Obtenida');


% Calcular y establecer los límites del gráfico
puntos_validos_plot = trayectoria_robot_final(~isnan(trayectoria_robot_final(:,1)), :);
x_min_plot=min(puntos_validos_plot(:,1)); x_max_plot=max(puntos_validos_plot(:,1));
y_min_plot=min(puntos_validos_plot(:,2)); y_max_plot=max(puntos_validos_plot(:,2));
padding = 2;
xlim([x_min_plot-padding, x_max_plot+padding]); ylim([y_min_plot-padding, y_max_plot+padding]);

% Crear objetos gráficos vacíos
h_brazo = plot(NaN, NaN, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Brazo');
h_real = plot(NaN, NaN, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Trayectoria Real');
legend('Location', 'northeastoutside');

% Encuentra el primer ángulo válido en tu historial
primer_angulo_idx = find(~cellfun(@(c) all(isnan(c)), historial_angulos), 1, 'first');
if isempty(primer_angulo_idx)
    error('No se encontraron ángulos válidos en el historial.');
end
theta_inicial = historial_angulos{primer_angulo_idx};

% --- PASO 2: Bucle de Dibujo (Comienza Directamente) ---

% Dibuja el brazo en su posición inicial (el primer punto de la figura)
[x0,y0,x1,y1,x2,y2,x3,y3] = calcular_pos_brazo(theta_inicial, L1, L2, L3);
set(h_brazo, 'XData', [x0, x1, x2, x3], 'YData', [y0, y1, y2, y3]);

% Inicializa la trayectoria real con un solo punto (la punta del brazo)
set(h_real, 'XData', x3, 'YData', y3);
drawnow;
pause(5.0); % Pausa para ver la posición inicial

% El bucle empieza desde el SIGUIENTE punto
for i = (primer_angulo_idx + 1):numel(historial_angulos)
    thetas = historial_angulos{i};
    
    if all(isnan(thetas))
        % Para los saltos, insertamos un NaN para levantar la pluma
        current_x = get(h_real, 'XData');
        current_y = get(h_real, 'YData');
        set(h_real, 'XData', [current_x, NaN], 'YData', [current_y, NaN]);
        continue;
    end
    
    % Mover el brazo
    [x0,y0,x1,y1,x2,y2,x3,y3] = calcular_pos_brazo(thetas, L1, L2, L3);
    set(h_brazo, 'XData', [x0, x1, x2, x3], 'YData', [y0, y1, y2, y3]);
    
    % Añadir la nueva posición de la punta a la trayectoria ya existente
    current_x = get(h_real, 'XData');
    current_y = get(h_real, 'YData');
    set(h_real, 'XData', [current_x, x3], 'YData', [current_y, y3]);
    
    drawnow;
    pause(0.01);
end

fprintf('✅ Animación completada.\n');

% (La función auxiliar 'calcular_pos_brazo' se mantiene igual)

% --- Función auxiliar ---
function [x0,y0,x1,y1,x2,y2,x3,y3] = calcular_pos_brazo(thetas, L1, L2, L3)
    th1 = thetas(1); th2 = thetas(2); th3 = thetas(3);
    x0 = 0; y0 = 0;
    x1 = L1*cos(th1); y1 = L1*sin(th1);
    x2 = x1 + L2*cos(th1+th2); y2 = y1 + L2*sin(th1+th2);
    x3 = x2 + L3*cos(th1+th2+th3); y3 = y2 + L3*sin(th1+th2+th3);
end
%% ============================================================
%  PARTE 4: PREPARACIÓN DE DATOS PARA ANÁLISIS NUMÉRICO
% =============================================================

% Vamos a trabajar con el primer trazo como ejemplo.
% 'historial_angulos' es una celda. Busquemos el primer segmento continuo.
idx_fin_primer_trazo = find(cellfun(@(c) all(isnan(c)), historial_angulos), 1, 'first') - 1;
if isempty(idx_fin_primer_trazo)
    idx_fin_primer_trazo = numel(historial_angulos); % Por si solo hay un trazo
end

% Extraer los datos del primer trazo
puntos_primer_trazo = trayectoria_dibujada(1:idx_fin_primer_trazo, :);
angulos_primer_trazo_cell = historial_angulos(1:idx_fin_primer_trazo);
angulos_primer_trazo = cell2mat(angulos_primer_trazo_cell'); % Convertir a matriz Nx3

% Crear un vector de tiempo simple. Asumamos que cada paso toma 0.1 segundos.
dt = 0.1; % Delta de tiempo entre puntos
num_puntos = size(puntos_primer_trazo, 1);
tiempo = (0:num_puntos-1)' * dt;

% Datos listos para la interpolación:
% tiempo -> vector columna de tiempo
% puntos_primer_trazo(:, 1) -> vector de X
% puntos_primer_trazo(:, 2) -> vector de Y
% angulos_primer_trazo(:, 1) -> vector de theta1
% etc.

% --- Función de Interpolación de Lagrange (implementación propia) ---
function yi = lagrange_interp(x, y, xi)
    n = length(x);
    yi = 0;
    for i = 1:n
        % Calcular el polinomio base L_i(xi)
        L = 1;
        for j = 1:n
            if i ~= j
                L = L * (xi - x(j)) / (x(i) - x(j));
            end
        end
        yi = yi + y(i) * L;
    end
end

%% ============================================================
%  PARTE 5: APLICACIÓN Y ANÁLISIS DE INTERPOLACIÓN (CON SPLINE)
% =============================================================

% (Se asume que ya tienes 'tiempo' y 'puntos_primer_trazo' de antes)

% Puntos de evaluación más finos para una gráfica suave
tiempo_fino = linspace(min(tiempo), max(tiempo), 500)';

% --- Aplicar Spline Cúbico ---
% 1. Construir el spline: calcular los coeficientes una sola vez
coeffs_spline_x = build_natural_cubic_spline(tiempo, puntos_primer_trazo(:, 1));
% 2. Evaluar el spline en los puntos finos
X_spline = eval_natural_cubic_spline(coeffs_spline_x, tiempo_fino);


% --- (Opcional) Re-calcular Lagrange para la comparación ---
% Si tienes muchos puntos, esto puede tardar mucho o fallar.
% X_lagrange = zeros(size(tiempo_fino));
% for k = 1:length(tiempo_fino)
%     X_lagrange(k) = lagrange_interp(tiempo, puntos_primer_trazo(:, 1), tiempo_fino(k));
% end

% --- Visualización Comparativa ---
figure;
hold on;
plot(tiempo, puntos_primer_trazo(:, 1), 'bo', 'DisplayName', 'Puntos Originales X(t)');
plot(tiempo_fino, X_spline, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Spline Cúbico');

% Descomenta la siguiente línea para ver el desastre de Lagrange de nuevo
% plot(tiempo_fino, X_lagrange, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Lagrange');

legend;
title('Comparación de Interpolación para X(t)');
xlabel('Tiempo (s)');
ylabel('Posición X (cm)');
grid on;

% Ajusta los límites del eje Y para ignorar los valores extremos de Lagrange
ylim([min(puntos_primer_trazo(:, 1)) - 1, max(puntos_primer_trazo(:, 1)) + 1]);
%%
% =========================================================================
%  FUNCIONES DE SPLINE CÚBICO NATURAL
% =========================================================================

function spline_coeffs = build_natural_cubic_spline(x, y)
    % Esta función calcula los coeficientes (segundas derivadas) para un
    % spline cúbico natural.
    % ENTRADA:
    %   x: vector de coordenadas x de los datos (debe estar ordenado)
    %   y: vector de coordenadas y de los datos
    % SALIDA:
    %   spline_coeffs: una estructura que contiene todos los coeficientes
    %                  necesarios para la evaluación.

    n = length(x);
    
    % Paso 1: Calcular los 'h', que son las distancias entre los puntos x
    h = zeros(n-1, 1);
    for i = 1:n-1
        h(i) = x(i+1) - x(i);
    end

    % Paso 2: Configurar el sistema de ecuaciones tridiagonal A*M = b
    % para encontrar las segundas derivadas M.
    % Para un spline NATURAL, M(1) = 0 y M(n) = 0, por lo que solo
    % resolvemos para los puntos interiores M(2)...M(n-1).
    
    % Crear la matriz A (tamaño n-2 x n-2)
    A = zeros(n-2, n-2);
    for i = 1:n-2
        A(i,i) = 2 * (h(i) + h(i+1));
        if i > 1
            A(i, i-1) = h(i);
        end
        if i < n-2
            A(i, i+1) = h(i+1);
        end
    end

    % Crear el vector b
    b = zeros(n-2, 1);
    for i = 1:n-2
        b(i) = 6 * ( (y(i+2)-y(i+1))/h(i+1) - (y(i+1)-y(i))/h(i) );
    end
    
    % Paso 3: Resolver el sistema para las segundas derivadas interiores
    M_internal = A\b; % El operador '\' de MATLAB es un solver eficiente
    
    % Paso 4: Construir el vector completo de segundas derivadas,
    % incluyendo los extremos (que son 0 para un spline natural).
    M = [0; M_internal; 0];
    
    % Guardar todo en una estructura para pasarlo a la función de evaluación
    spline_coeffs.x = x;
    spline_coeffs.y = y;
    spline_coeffs.M = M;
    spline_coeffs.h = h;
end


function y_eval = eval_natural_cubic_spline(spline_coeffs, x_eval)
    % Esta función evalúa el spline en un nuevo conjunto de puntos x_eval.
    % ENTRADA:
    %   spline_coeffs: la estructura generada por build_natural_cubic_spline
    %   x_eval: los puntos donde se quiere evaluar el spline
    % SALIDA:
    %   y_eval: los valores interpolados del spline
    
    % Extraer los coeficientes de la estructura
    x = spline_coeffs.x;
    y = spline_coeffs.y;
    M = spline_coeffs.M;
    h = spline_coeffs.h;
    
    n = length(x);
    y_eval = zeros(size(x_eval));
    
    % Para cada punto a evaluar...
    for k = 1:length(x_eval)
        xk = x_eval(k);
        
        % Paso 1: Encontrar en qué intervalo [x_i, x_{i+1}] se encuentra xk
        % Buscamos el último x que es menor o igual a xk
        i = find(x <= xk, 1, 'last');
        
        % Manejar el caso del último punto para evitar errores de índice
        if i == n
            i = n - 1;
        end
        
        % Paso 2: Aplicar la fórmula de evaluación del spline para ese intervalo
        % Esta es la fórmula estándar del polinomio cúbico del spline
        A = (M(i+1) - M(i)) / (6*h(i));
        B = M(i) / 2;
        C = (y(i+1) - y(i))/h(i) - h(i)*(M(i+1) + 2*M(i))/6;
        D = y(i);
        
        t = xk - x(i); % Término de distancia
        
        yk = A*t^3 + B*t^2 + C*t + D; % Polinomio de Horner sería más eficiente, pero esto es más claro
        
        y_eval(k) = yk;
    end
end