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

% Vamos a trabajar con el primer trazo continuo como ejemplo.
% 'historial_angulos' es una celda. Buscamos dónde termina el primer segmento.
idx_fin_primer_trazo = find(cellfun(@(c) all(isnan(c)), historial_angulos), 1, 'first') - 1;
if isempty(idx_fin_primer_trazo)
    idx_fin_primer_trazo = numel(historial_angulos); % En caso de que solo haya un trazo.
end

% --- Extraer los datos del primer trazo ---
% Puntos objetivo (X, Y) del primer trazo.
puntos_primer_trazo = trayectoria_robot_final(1:idx_fin_primer_trazo, :);

% Ángulos calculados por la cinemática inversa para ese trazo.
angulos_primer_trazo_cell = historial_angulos(1:idx_fin_primer_trazo);

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% CORRECCIÓN CLAVE AQUÍ: Transponer DESPUÉS de cell2mat
% Esto convierte la celda de vectores columna [3x1] en una matriz Nx3.
angulos_primer_trazo = cell2mat(angulos_primer_trazo_cell)'; 
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Crear un vector de tiempo simple para usar como variable independiente.
dt = 0.1; % Delta de tiempo entre puntos (ej. 0.1 segundos)
num_puntos = size(puntos_primer_trazo, 1);
tiempo = (0:num_puntos-1)' * dt;

% Puntos de evaluación más finos para una gráfica suave en las interpolaciones
tiempo_fino = linspace(min(tiempo), max(tiempo), 500)';

fprintf('✅ Datos del primer trazo preparados para el análisis numérico.\n');
fprintf('   - Número de puntos en el trazo: %d\n', num_puntos);
fprintf('   - Duración simulada del trazo: %.2f segundos\n', max(tiempo));


%% ============================================================
%  PARTE 5: ANÁLISIS CON INTERPOLACIÓN DE LAGRANGE
% =============================================================
% ADVERTENCIA: La interpolación de Lagrange con muchos puntos es numéricamente
% inestable y sufre del fenómeno de Runge (oscilaciones salvajes en los bordes).
% Lo implementamos aquí con fines educativos y de comparación.

fprintf('\nIniciando interpolación con Lagrange (esto puede ser lento)... \n');

% Pre-alocar memoria para los resultados
X_lagrange = zeros(size(tiempo_fino));
Y_lagrange = zeros(size(tiempo_fino));

% Evaluar la interpolación en cada punto de 'tiempo_fino'
for k = 1:length(tiempo_fino)
    X_lagrange(k) = lagrange_interp(tiempo, puntos_primer_trazo(:, 1), tiempo_fino(k));
    Y_lagrange(k) = lagrange_interp(tiempo, puntos_primer_trazo(:, 2), tiempo_fino(k));
end

fprintf('✅ Interpolación con Lagrange completada.\n');

% --- Visualización de los resultados de Lagrange ---
figure('Name', 'Análisis de Interpolación de Lagrange');
plot(puntos_primer_trazo(:, 1), puntos_primer_trazo(:, 2), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Puntos Originales');
hold on;
plot(X_lagrange, Y_lagrange, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Interpolación Lagrange');
grid on;
axis equal;
title('Resultado de la Interpolación de Lagrange');
xlabel('Posición X (cm)');
ylabel('Posición Y (cm)');
legend;
padding = 2;
xlim([min(puntos_primer_trazo(:, 1)) - padding, max(puntos_primer_trazo(:, 1)) + padding]);
ylim([min(puntos_primer_trazo(:, 2)) - padding, max(puntos_primer_trazo(:, 2)) + padding]);


%% ============================================================
%  PARTE 6: ANÁLISIS CON INTERPOLACIÓN DE SPLINES CÚBICOS
% =============================================================
% Los splines cúbicos son mucho más estables y adecuados para trayectorias.

fprintf('\nIniciando interpolación con Splines Cúbicos...\n');

% --- Interpolar la trayectoria (X, Y) ---
coeffs_spline_x = build_natural_cubic_spline(tiempo, puntos_primer_trazo(:, 1));
coeffs_spline_y = build_natural_cubic_spline(tiempo, puntos_primer_trazo(:, 2));

X_spline = eval_natural_cubic_spline(coeffs_spline_x, tiempo_fino);
Y_spline = eval_natural_cubic_spline(coeffs_spline_y, tiempo_fino);


% --- Interpolar los ángulos (Theta1, Theta2, Theta3) ---
coeffs_th1 = build_natural_cubic_spline(tiempo, angulos_primer_trazo(:, 1));
coeffs_th2 = build_natural_cubic_spline(tiempo, angulos_primer_trazo(:, 2));
coeffs_th3 = build_natural_cubic_spline(tiempo, angulos_primer_trazo(:, 3));

theta1_spline = eval_natural_cubic_spline(coeffs_th1, tiempo_fino);
theta2_spline = eval_natural_cubic_spline(coeffs_th2, tiempo_fino);
theta3_spline = eval_natural_cubic_spline(coeffs_th3, tiempo_fino);

fprintf('✅ Interpolación con Splines completada.\n');


% --- Visualización Comparativa (Lagrange vs. Spline) ---
figure('Name', 'Comparación de Métodos de Interpolación');
subplot(2, 2, 1);
hold on;
plot(tiempo, puntos_primer_trazo(:, 1), 'ko', 'DisplayName', 'Puntos Originales');
plot(tiempo_fino, X_spline, 'g-', 'LineWidth', 2, 'DisplayName', 'Spline Cúbico');
plot(tiempo_fino, X_lagrange, 'r--', 'DisplayName', 'Lagrange');
title('Comparación para X(t)'); xlabel('Tiempo (s)'); ylabel('Posición X (cm)');
legend('Location', 'best'); grid on;
ylim_padding = 1;
ylim([min(puntos_primer_trazo(:, 1)) - ylim_padding, max(puntos_primer_trazo(:, 1)) + ylim_padding]);

subplot(2, 2, 2);
hold on;
plot(tiempo, puntos_primer_trazo(:, 2), 'ko', 'DisplayName', 'Puntos Originales');
plot(tiempo_fino, Y_spline, 'g-', 'LineWidth', 2, 'DisplayName', 'Spline Cúbico');
plot(tiempo_fino, Y_lagrange, 'r--', 'DisplayName', 'Lagrange');
title('Comparación para Y(t)'); xlabel('Tiempo (s)'); ylabel('Posición Y (cm)');
legend('Location', 'best'); grid on;
ylim([min(puntos_primer_trazo(:, 2)) - ylim_padding, max(puntos_primer_trazo(:, 2)) + ylim_padding]);

subplot(2, 2, 3);
hold on;
plot(tiempo, rad2deg(angulos_primer_trazo(:, 1)), 'bo', 'DisplayName', 'Puntos Originales');
plot(tiempo_fino, rad2deg(theta1_spline), 'c-', 'LineWidth', 2, 'DisplayName', 'Spline Cúbico');
title('Interpolación de Ángulo \theta_1(t)'); xlabel('Tiempo (s)'); ylabel('Ángulo (grados)');
legend('Location', 'best'); grid on;

subplot(2, 2, 4);
hold on;
plot(X_spline, Y_spline, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Trayectoria Spline');
plot(puntos_primer_trazo(:, 1), puntos_primer_trazo(:, 2), 'k.', 'MarkerSize', 10, 'DisplayName', 'Puntos Originales');
axis equal; grid on;
title('Trayectoria 2D Reconstruida con Splines'); xlabel('Posición X (cm)'); ylabel('Posición Y (cm)');
legend('Location', 'best');

%% =========================================================================
%  PARTE 7: FUNCIONES AUXILIARES DE INTERPOLACIÓN
% =========================================================================
% (Esta parte no necesita cambios, la incluyo para que tengas el bloque completo)

% --- Función de Interpolación de Lagrange ---
function yi = lagrange_interp(x, y, xi)
    n = length(x);
    yi = 0;
    for i = 1:n
        L = 1;
        for j = 1:n
            if i ~= j
                L = L * (xi - x(j)) / (x(i) - x(j));
            end
        end
        yi = yi + y(i) * L;
    end
end

% --- Funciones de Spline Cúbico Natural ---
function spline_coeffs = build_natural_cubic_spline(x, y)
    n = length(x);
    if n < 3
        error('Se necesitan al menos 3 puntos para un spline cúbico.');
    end
    h = zeros(n-1, 1);
    for i = 1:n-1
        h(i) = x(i+1) - x(i);
    end
    A = zeros(n-2, n-2);
    b = zeros(n-2, 1);
    for i = 1:n-2
        A(i,i) = 2 * (h(i) + h(i+1));
        if i > 1, A(i, i-1) = h(i); end
        if i < n-2, A(i, i+1) = h(i+1); end
        term1 = (y(i+2) - y(i+1)) / h(i+1);
        term2 = (y(i+1) - y(i)) / h(i);
        b(i) = 6 * (term1 - term2);
    end
    M_internal = A\b;
    M = [0; M_internal; 0];
    spline_coeffs.x = x;
    spline_coeffs.y = y;
    spline_coeffs.M = M;
    spline_coeffs.h = h;
end

function y_eval = eval_natural_cubic_spline(spline_coeffs, x_eval)
    x = spline_coeffs.x; y = spline_coeffs.y;
    M = spline_coeffs.M; h = spline_coeffs.h;
    n = length(x);
    y_eval = zeros(size(x_eval));
    for k = 1:length(x_eval)
        xk = x_eval(k);
        i = find(x <= xk, 1, 'last');
        if i == n, i = n - 1; end
        t = xk - x(i);
        A = (M(i+1) - M(i)) / (6*h(i));
        B = M(i) / 2;
        C = (y(i+1) - y(i))/h(i) - h(i)*(2*M(i) + M(i+1))/6;
        D = y(i);
        y_eval(k) = A*t^3 + B*t^2 + C*t + D;
    end
end

%% =========================================================================
%  PARTE 8: DIFERENCIACIÓN NUMÉRICA (CÁLCULO DE VELOCIDAD Y ACELERACIÓN)
% =========================================================================
% Objetivo: A partir de la trayectoria suave obtenida con splines (posición vs. tiempo),
% vamos a calcular la velocidad y la aceleración en cada punto.
% Usaremos el método de diferencias finitas centrales, que es más preciso.
% v(t) ≈ (p(t+h) - p(t-h)) / (2h)
% a(t) ≈ (p(t+h) - 2p(t) + p(t-h)) / h^2

fprintf('\nIniciando Parte 3: Diferenciación e Integración Numérica...\n');

% --- Calcular velocidad y aceleración para la trayectoria cartesiana (X, Y) ---
[velocidad_x, aceleracion_x] = diferenciacion_numerica(tiempo_fino, X_spline);
[velocidad_y, aceleracion_y] = diferenciacion_numerica(tiempo_fino, Y_spline);

% --- Calcular velocidad y aceleración para los ángulos de las juntas (Thetas) ---
[vel_angular_1, acel_angular_1] = diferenciacion_numerica(tiempo_fino, theta1_spline);
[vel_angular_2, acel_angular_2] = diferenciacion_numerica(tiempo_fino, theta2_spline);
[vel_angular_3, acel_angular_3] = diferenciacion_numerica(tiempo_fino, theta3_spline);

fprintf('✅ Diferenciación numérica completada. Velocidades y aceleraciones calculadas.\n');

% --- Visualización de las derivadas ---
figure('Name', 'Análisis Cinemático (Velocidad y Aceleración)');
subplot(2, 2, 1);
plot(tiempo_fino, velocidad_x, 'b-', 'DisplayName', 'Velocidad X');
hold on;
plot(tiempo_fino, velocidad_y, 'r-', 'DisplayName', 'Velocidad Y');
title('Velocidad Cartesiana vs. Tiempo');
xlabel('Tiempo (s)'); ylabel('Velocidad (cm/s)'); grid on; legend;

subplot(2, 2, 2);
plot(tiempo_fino, aceleracion_x, 'b-', 'DisplayName', 'Aceleración X');
hold on;
plot(tiempo_fino, aceleracion_y, 'r-', 'DisplayName', 'Aceleración Y');
title('Aceleración Cartesiana vs. Tiempo');
xlabel('Tiempo (s)'); ylabel('Aceleración (cm/s^2)'); grid on; legend;

subplot(2, 2, 3);
plot(tiempo_fino, rad2deg(vel_angular_1), 'DisplayName', '\omega_1');
hold on;
plot(tiempo_fino, rad2deg(vel_angular_2), 'DisplayName', '\omega_2');
plot(tiempo_fino, rad2deg(vel_angular_3), 'DisplayName', '\omega_3');
title('Velocidad Angular vs. Tiempo');
xlabel('Tiempo (s)'); ylabel('Velocidad Angular (deg/s)'); grid on; legend;

subplot(2, 2, 4);
plot(tiempo_fino, rad2deg(acel_angular_1), 'DisplayName', '\alpha_1');
hold on;
plot(tiempo_fino, rad2deg(acel_angular_2), 'DisplayName', '\alpha_2');
plot(tiempo_fino, rad2deg(acel_angular_3), 'DisplayName', '\alpha_3');
title('Aceleración Angular vs. Tiempo');
xlabel('Tiempo (s)'); ylabel('Aceleración Angular (deg/s^2)'); grid on; legend;


%% =========================================================================
%  PARTE 9: INTEGRACIÓN NUMÉRICA (CÁLCULO DE IMPULSO Y TRABAJO)
% =========================================================================
% Objetivo: Aplicar la Regla de Simpson y la Cuadratura Gaussiana para
% calcular el impulso total y el trabajo realizado por los motores.

% --- Definición de Parámetros Físicos (según el enunciado y supuestos) ---
% Convertimos todo a unidades del Sistema Internacional (metros, kg, segundos)
masa_marcador = 3 / 1000; % 3g a kg
masa_soporte = 10 / 1000; % 10g a kg
masa_efector_final = masa_marcador + masa_soporte; % Masa total en la punta

% Para el torque, necesitamos la inercia de los eslabones. Como no se da,
% la modelamos como una varilla delgada que rota por un extremo: I = (1/3)mL^2
% Asumimos una masa para los eslabones (ej. 100g para los largos, 50g para el corto)
masa_L1 = 0.1; % kg
masa_L2 = 0.1; % kg
masa_L3 = 0.05; % kg

L1_m = L1 / 100; % cm a m
L2_m = L2 / 100; % cm a m
L3_m = L3 / 100; % cm a m

% Momentos de inercia aproximados
I1 = (1/3) * masa_L1 * L1_m^2;
I2 = (1/3) * masa_L2 * L2_m^2;
I3 = (1/3) * masa_L3 * L3_m^2;

% --- Cálculo de Fuerza y Torque ---
% Fuerza F = m*a (en Newtons)
fuerza_x = masa_efector_final * (aceleracion_x / 100); % Aceleración de cm/s^2 a m/s^2
fuerza_y = masa_efector_final * (aceleracion_y / 100);

% Torque T = I*alpha (en Newton-metro)
torque_1 = I1 * acel_angular_1;
torque_2 = I2 * acel_angular_2;
torque_3 = I3 * acel_angular_3;


fprintf('\nCalculando integrales numéricas...\n');
% --- Cálculo del Impulso Total (Integral de Fuerza vs Tiempo) ---
% Impulso en X
impulso_x_simpson = regla_simpson(tiempo_fino, fuerza_x);
impulso_x_gauss = cuadratura_gaussiana(tiempo_fino, fuerza_x);
% Impulso en Y
impulso_y_simpson = regla_simpson(tiempo_fino, fuerza_y);
impulso_y_gauss = cuadratura_gaussiana(tiempo_fino, fuerza_y);

% --- Cálculo del Trabajo Total (Integral de Torque vs Angulo) ---
% Trabajo del motor 1
trabajo_1_simpson = regla_simpson(theta1_spline, torque_1);
trabajo_1_gauss = cuadratura_gaussiana(theta1_spline, torque_1);
% Trabajo del motor 2
trabajo_2_simpson = regla_simpson(theta2_spline, torque_2);
trabajo_2_gauss = cuadratura_gaussiana(theta2_spline, torque_2);
% Trabajo del motor 3
trabajo_3_simpson = regla_simpson(theta3_spline, torque_3);
trabajo_3_gauss = cuadratura_gaussiana(theta3_spline, torque_3);


% --- Mostrar Resultados en Consola ---
fprintf('----------------- RESULTADOS DE INTEGRACIÓN -----------------\n');
fprintf('Impulso Total en X (Simpson): %.6f N·s\n', impulso_x_simpson);
fprintf('Impulso Total en X (Gauss):   %.6f N·s\n', impulso_x_gauss);
fprintf('Impulso Total en Y (Simpson): %.6f N·s\n', impulso_y_simpson);
fprintf('Impulso Total en Y (Gauss):   %.6f N·s\n', impulso_y_gauss);
fprintf('------------------------------------------------------------\n');
fprintf('Trabajo Motor 1 (Simpson): %.6f Joules\n', trabajo_1_simpson);
fprintf('Trabajo Motor 1 (Gauss):   %.6f Joules\n', trabajo_1_gauss);
fprintf('Trabajo Motor 2 (Simpson): %.6f Joules\n', trabajo_2_simpson);
fprintf('Trabajo Motor 2 (Gauss):   %.6f Joules\n', trabajo_2_gauss);
fprintf('Trabajo Motor 3 (Simpson): %.6f Joules\n', trabajo_3_simpson);
fprintf('Trabajo Motor 3 (Gauss):   %.6f Joules\n', trabajo_3_gauss);
fprintf('------------------------------------------------------------\n');


%% =========================================================================
%  PARTE 10: VISUALIZACIÓN DE RESULTADOS DE INTEGRACIÓN
% =========================================================================
% Objetivo: Graficar la magnitud del impulso acumulado a lo largo del tiempo.

% Para esto, calculamos la integral acumulativa. La función 'cumtrapz' de
% MATLAB es útil para esto (integral por trapecio acumulada).
impulso_acum_x = cumtrapz(tiempo_fino, fuerza_x);
impulso_acum_y = cumtrapz(tiempo_fino, fuerza_y);

% Magnitud del vector de impulso en cada instante de tiempo
magnitud_impulso_acum = sqrt(impulso_acum_x.^2 + impulso_acum_y.^2);

figure('Name', 'Análisis de Impulso');
plot(tiempo_fino, magnitud_impulso_acum, 'm-', 'LineWidth', 2);
title('Magnitud del Impulso Acumulado vs. Tiempo');
xlabel('Tiempo (s)');
ylabel('Magnitud del Impulso |J| (N·s)');
grid on;
fprintf('✅ Visualización del impulso completada.\n');


%% =========================================================================
%  PARTE 11: FUNCIONES AUXILIARES PARA DIFERENCIACIÓN E INTEGRACIÓN
% =========================================================================
% (Añadir estas funciones al final de tu script, junto con las demás)

function [v, a] = diferenciacion_numerica(t, p)
    % Calcula la primera y segunda derivada de p(t) usando diferencias centrales.
    n = length(t);
    h = t(2) - t(1); % Asumimos paso de tiempo constante
    v = zeros(n, 1);
    a = zeros(n, 1);

    % Puntos interiores (más precisos)
    for i = 2:n-1
        v(i) = (p(i+1) - p(i-1)) / (2*h);
        a(i) = (p(i+1) - 2*p(i) + p(i-1)) / (h^2);
    end

    % Extremos (menos precisos, usando diferencias hacia adelante/atrás)
    v(1) = (p(2) - p(1)) / h;
    v(n) = (p(n) - p(n-1)) / h;
    a(1) = (v(2) - v(1)) / h;
    a(n) = (v(n) - v(n-1)) / h;
end

function integral = regla_simpson(x, y)
    % Implementa la Regla Compuesta de Simpson 1/3.
    % x: vector de puntos (variable independiente), debe ser equidistante.
    % y: vector de valores de la función a integrar.
    n = length(x);
    if mod(n-1, 2) ~= 0
        warning('El número de intervalos para Simpson no es par. El último intervalo se calculará con la regla del trapecio.');
        % En un caso real, se podría interpolar un punto extra o truncar.
        % Aquí, para robustez, calcularemos Simpson hasta n-1 y añadimos el último trapecio.
        integral_simpson = regla_simpson(x(1:n-1), y(1:n-1));
        integral_trapecio = (x(n)-x(n-1)) * (y(n)+y(n-1))/2;
        integral = integral_simpson + integral_trapecio;
        return;
    end
    
    h = x(2) - x(1); % Asumimos paso constante
    
    suma_impares = sum(y(2:2:n-1));
    suma_pares = sum(y(3:2:n-2));
    
    integral = (h/3) * (y(1) + 4*suma_impares + 2*suma_pares + y(n));
end

function integral = cuadratura_gaussiana(x_data, y_data)
    % Implementa la Cuadratura de Gauss-Legendre de 2 puntos sobre intervalos.
    % x_data, y_data: Puntos de datos. No necesitan ser equidistantes.
    
    integral = 0;
    n = length(x_data);
    
    % Puntos y pesos para Gauss-Legendre de 2 puntos
    w1 = 1.0; w2 = 1.0;
    t1 = -1/sqrt(3); t2 = 1/sqrt(3);

    % Iteramos sobre cada intervalo [x_i, x_{i+1}]
    for i = 1:n-1
        a = x_data(i);
        b = x_data(i+1);
        
        % Transformación de la variable de integración de [a,b] a [-1,1]
        % x = (b-a)/2 * t + (a+b)/2
        % dx = (b-a)/2 * dt
        
        % Puntos de evaluación en el intervalo original
        x1 = (b-a)/2 * t1 + (a+b)/2;
        x2 = (b-a)/2 * t2 + (a+b)/2;
        
        % Evaluamos la función en esos puntos. Como no tenemos la función analítica,
        % usamos una interpolación lineal simple entre los puntos de datos.
        y1 = interp1(x_data, y_data, x1);
        y2 = interp1(x_data, y_data, x2);
        
        % Sumamos la contribución del intervalo a la integral total
        integral = integral + (b-a)/2 * (w1*y1 + w2*y2);
    end
end


%% =========================================================================
%  PARTE 12: BONO - CINEMÁTICA INVERSA CON OPTIMIZACIÓN (NEWTON-RAPHSON 3x3)
% =========================================================================
% Objetivo: Resolver el sistema cinemático añadiendo una tercera ecuación de
% restricción para hacerlo un sistema cuadrado 3x3 y resolverlo con el
% método de Newton-Raphson estándar.
%
% Restricción elegida: Mantener la orientación del efector final constante.
% Ecuación: theta1 + theta2 + theta3 = 0 (último eslabón horizontal).

fprintf('\nIniciando Parte de Bono: Solver con Newton-Raphson 3x3...\n');

% Posición inicial del robot (la misma que antes)
theta_inicial_deg = [10; 0; 0]; 
theta_actual_rad_NR = deg2rad(theta_inicial_deg);

% Almacenaremos los nuevos ángulos calculados aquí
historial_angulos_optimizado = {};

% Recorrer cada trazo de la trayectoria
for i = 1:size(trayectoria_robot_final, 1)
    
    punto_actual = trayectoria_robot_final(i, :);
    
    % Si es un separador de trazos (NaN), lo guardamos y continuamos
    if isnan(punto_actual(1))
        historial_angulos_optimizado{end+1} = NaN(3,1);
        continue;
    end
    
    % Punto objetivo actual [x; y]
    punto_objetivo = [punto_actual(1); punto_actual(2)];
    
    % Resolver la cinemática inversa con el nuevo solver de Newton-Raphson
    theta_actual_rad_NR = inverse_kinematics_solver_NR(punto_objetivo, theta_actual_rad_NR, L1, L2, L3);
    
    % Guardar el resultado
    historial_angulos_optimizado{end+1} = theta_actual_rad_NR;
end

fprintf('✅ ¡Cálculo de cinemática inversa con optimización completado!\n');


%% =========================================================================
%  PARTE 13: COMPARACIÓN VISUAL DE AMBOS MÉTODOS
% =========================================================================

figure('Name', 'Comparación de Soluciones Cinemáticas');
hold on; axis equal; grid on;
title('Solución Original (Pseudoinversa) vs. Solución Optimizada (Newton-Raphson 3x3)');
xlabel('X (cm)'); ylabel('Y (cm)');

% Dibujar la trayectoria objetivo
plot(trayectoria_robot_final(:,1), trayectoria_robot_final(:,2), 'k:', 'LineWidth', 2, 'DisplayName', 'Trayectoria Objetivo');

% Encontrar el primer ángulo válido en ambos historiales
idx_start_orig = find(~cellfun(@(c) all(isnan(c)), historial_angulos), 1, 'first');
idx_start_opt = find(~cellfun(@(c) all(isnan(c)), historial_angulos_optimizado), 1, 'first');

% Calcular la posición inicial de ambos brazos
[~,~,~,~,~,~,x3_orig,y3_orig] = calcular_pos_brazo(historial_angulos{idx_start_orig}, L1, L2, L3);
[~,~,~,~,~,~,x3_opt,y3_opt] = calcular_pos_brazo(historial_angulos_optimizado{idx_start_opt}, L1, L2, L3);

% Crear objetos gráficos para las trayectorias dibujadas
h_real_orig = plot(x3_orig, y3_orig, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Dibujo Original (pinv)');
h_real_opt = plot(x3_opt, y3_opt, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Dibujo Optimizado (NR 3x3)');
legend('Location', 'northeastoutside');

% Extraer todas las posiciones de la punta del efector para ambos métodos
puntos_dibujados_orig = cellfun(@(th) forward_kinematics(th, L1, L2, L3), historial_angulos, 'UniformOutput', false);
puntos_dibujados_orig = cell2mat(puntos_dibujados_orig)';

puntos_dibujados_opt = cellfun(@(th) forward_kinematics(th, L1, L2, L3), historial_angulos_optimizado, 'UniformOutput', false);
puntos_dibujados_opt = cell2mat(puntos_dibujados_opt)';

% Actualizar las gráficas con las trayectorias completas
set(h_real_orig, 'XData', puntos_dibujados_orig(:,1), 'YData', puntos_dibujados_orig(:,2));
set(h_real_opt, 'XData', puntos_dibujados_opt(:,1), 'YData', puntos_dibujados_opt(:,2));

% --- Gráfica de comparación de ángulos ---
angulos_orig = cell2mat(historial_angulos)';
angulos_opt = cell2mat(historial_angulos_optimizado)';

figure('Name', 'Comparación de Ángulos');
subplot(3,1,1);
plot(rad2deg(angulos_orig(:,1)), 'b-'); hold on;
plot(rad2deg(angulos_opt(:,1)), 'r-');
title('\theta_1'); ylabel('Grados'); grid on; legend('Original', 'Optimizado');

subplot(3,1,2);
plot(rad2deg(angulos_orig(:,2)), 'b-'); hold on;
plot(rad2deg(angulos_opt(:,2)), 'r-');
title('\theta_2'); ylabel('Grados'); grid on;

subplot(3,1,3);
plot(rad2deg(angulos_orig(:,3)), 'b-'); hold on;
plot(rad2deg(angulos_opt(:,3)), 'r-');
title('\theta_3'); ylabel('Grados'); grid on;
xlabel('Índice del Punto');


%% =========================================================================
%  PARTE 14: FUNCIONES AUXILIARES PARA EL BONO
% =========================================================================
% (Añadir estas funciones al final de tu script, junto con las demás)

function J = jacobian_matrix_3x3(thetas, L1, L2, L3)
    th1 = thetas(1); th2 = thetas(2); th3 = thetas(3);
    J11 = -L1*sin(th1) - L2*sin(th1+th2) - L3*sin(th1+th2+th3);
    J12 = -L2*sin(th1+th2) - L3*sin(th1+th2+th3);
    J13 = -L3*sin(th1+th2+th3);
    J21 = L1*cos(th1) + L2*cos(th1+th2) + L3*cos(th1+th2+th3);
    J22 = L2*cos(th1+th2) + L3*cos(th1+th2+th3);
    J23 = L3*cos(th1+th2+th3);
    J = [J11, J12, J13; J21, J22, J23; 1, 1, 1];
end


function thetas_final = inverse_kinematics_solver_NR(target_pos, current_thetas, L1, L2, L3)
    % Solver de IK usando Newton-Raphson con Damped Least Squares (DLS)
    thetas = current_thetas;
    tolerancia = 0.01; max_iteraciones = 100; lambda = 0.1;
    target_full = [target_pos; 0]; 
    for i = 1:max_iteraciones
        current_pos = forward_kinematics(thetas, L1, L2, L3);
        current_constraint = thetas(1) + thetas(2) + thetas(3);
        error_vector = target_full - [current_pos; current_constraint];
        if norm(error_vector) < tolerancia, break; end
        J_F = jacobian_matrix_3x3(thetas, L1, L2, L3);
        JtJ = J_F' * J_F;
        damping_matrix = lambda^2 * eye(3);
        delta_thetas = (JtJ + damping_matrix) \ (J_F' * error_vector);
        thetas = thetas + delta_thetas;
    end
    thetas_final = thetas;
end

%% =========================================================================
%  BONO (0.2): ANÁLISIS REPETIDO PARA LA SOLUCIÓN OPTIMIZADA
% =========================================================================
% Objetivo: Repetir las partes 2 y 3 (interpolación, diferenciación e
% integración) utilizando los datos generados por el solver optimizado de
% Newton-Raphson 3x3 para comparar los resultados.

fprintf('\n\n--- INICIANDO ANÁLISIS REPETIDO PARA LA SOLUCIÓN OPTIMIZADA ---\n');

% --- PASO 1: Preparar los datos de la solución optimizada ---
% Vamos a trabajar con el primer trazo continuo de la solución optimizada.
% Es probable que sea más corto o diferente al original.
idx_fin_primer_trazo_opt = find(cellfun(@(c) all(isnan(c)), historial_angulos_optimizado), 1, 'first') - 1;
if isempty(idx_fin_primer_trazo_opt)
    idx_fin_primer_trazo_opt = numel(historial_angulos_optimizado);
end

% Extraer los datos de ese primer trazo
puntos_primer_trazo_opt = trayectoria_robot_final(1:idx_fin_primer_trazo_opt, :);
angulos_primer_trazo_cell_opt = historial_angulos_optimizado(1:idx_fin_primer_trazo_opt);
angulos_primer_trazo_opt = cell2mat(angulos_primer_trazo_cell_opt)';

% Crear un nuevo vector de tiempo para este segmento
dt = 0.1; % Mantenemos el mismo delta de tiempo
num_puntos_opt = size(puntos_primer_trazo_opt, 1);
tiempo_opt = (0:num_puntos_opt-1)' * dt;
tiempo_fino_opt = linspace(min(tiempo_opt), max(tiempo_opt), 500)';

fprintf('✅ Datos de la solución optimizada preparados para análisis.\n');

% --- PASO 2: Interpolación con Splines ---
fprintf('Iniciando interpolación con Splines para la solución optimizada...\n');

% Interpolar trayectoria (X, Y)
coeffs_spline_x_opt = build_natural_cubic_spline(tiempo_opt, puntos_primer_trazo_opt(:, 1));
coeffs_spline_y_opt = build_natural_cubic_spline(tiempo_opt, puntos_primer_trazo_opt(:, 2));
X_spline_opt = eval_natural_cubic_spline(coeffs_spline_x_opt, tiempo_fino_opt);
Y_spline_opt = eval_natural_cubic_spline(coeffs_spline_y_opt, tiempo_fino_opt);

% Interpolar ángulos (Thetas)
coeffs_th1_opt = build_natural_cubic_spline(tiempo_opt, angulos_primer_trazo_opt(:, 1));
coeffs_th2_opt = build_natural_cubic_spline(tiempo_opt, angulos_primer_trazo_opt(:, 2));
coeffs_th3_opt = build_natural_cubic_spline(tiempo_opt, angulos_primer_trazo_opt(:, 3));
theta1_spline_opt = eval_natural_cubic_spline(coeffs_th1_opt, tiempo_fino_opt);
theta2_spline_opt = eval_natural_cubic_spline(coeffs_th2_opt, tiempo_fino_opt);
theta3_spline_opt = eval_natural_cubic_spline(coeffs_th3_opt, tiempo_fino_opt);

% --- PASO 3: Diferenciación Numérica ---
fprintf('Calculando velocidades y aceleraciones para la solución optimizada...\n');

% Calcular velocidad y aceleración cartesiana
[velocidad_x_opt, aceleracion_x_opt] = diferenciacion_numerica(tiempo_fino_opt, X_spline_opt);
[velocidad_y_opt, aceleracion_y_opt] = diferenciacion_numerica(tiempo_fino_opt, Y_spline_opt);

% Calcular velocidad y aceleración angular
[vel_angular_1_opt, acel_angular_1_opt] = diferenciacion_numerica(tiempo_fino_opt, theta1_spline_opt);
[vel_angular_2_opt, acel_angular_2_opt] = diferenciacion_numerica(tiempo_fino_opt, theta2_spline_opt);
[vel_angular_3_opt, acel_angular_3_opt] = diferenciacion_numerica(tiempo_fino_opt, theta3_spline_opt);

% --- PASO 4: Integración Numérica ---
fprintf('Calculando impulso y trabajo para la solución optimizada...\n');

% Calcular Fuerza y Torque (los parámetros físicos son los mismos)
fuerza_x_opt = masa_efector_final * (aceleracion_x_opt / 100);
fuerza_y_opt = masa_efector_final * (aceleracion_y_opt / 100);
torque_1_opt = I1 * acel_angular_1_opt;
torque_2_opt = I2 * acel_angular_2_opt;
torque_3_opt = I3 * acel_angular_3_opt;

% Calcular Impulso Total
impulso_x_simpson_opt = regla_simpson(tiempo_fino_opt, fuerza_x_opt);
impulso_x_gauss_opt = cuadratura_gaussiana(tiempo_fino_opt, fuerza_x_opt);
impulso_y_simpson_opt = regla_simpson(tiempo_fino_opt, fuerza_y_opt);
impulso_y_gauss_opt = cuadratura_gaussiana(tiempo_fino_opt, fuerza_y_opt);

% Calcular Trabajo Total
trabajo_1_simpson_opt = regla_simpson(theta1_spline_opt, torque_1_opt);
trabajo_1_gauss_opt = cuadratura_gaussiana(theta1_spline_opt, torque_1_opt);
trabajo_2_simpson_opt = regla_simpson(theta2_spline_opt, torque_2_opt);
trabajo_2_gauss_opt = cuadratura_gaussiana(theta2_spline_opt, torque_2_opt);
trabajo_3_simpson_opt = regla_simpson(theta3_spline_opt, torque_3_opt);
trabajo_3_gauss_opt = cuadratura_gaussiana(theta3_spline_opt, torque_3_opt);

% --- PASO 5: Mostrar Resultados Numéricos y Gráficos ---

% Imprimir resultados de integración en consola
fprintf('\n---------- RESULTADOS DE INTEGRACIÓN (SOLUCIÓN OPTIMIZADA) ----------\n');
fprintf('Impulso Total en X (Simpson): %.6f N·s\n', impulso_x_simpson_opt);
fprintf('Impulso Total en X (Gauss):   %.6f N·s\n', impulso_x_gauss_opt);
fprintf('Impulso Total en Y (Simpson): %.6f N·s\n', impulso_y_simpson_opt);
fprintf('Impulso Total en Y (Gauss):   %.6f N·s\n', impulso_y_gauss_opt);
fprintf('---------------------------------------------------------------------\n');
fprintf('Trabajo Motor 1 (Simpson): %.6f Joules\n', trabajo_1_simpson_opt);
fprintf('Trabajo Motor 1 (Gauss):   %.6f Joules\n', trabajo_1_gauss_opt);
fprintf('Trabajo Motor 2 (Simpson): %.6f Joules\n', trabajo_2_simpson_opt);
fprintf('Trabajo Motor 2 (Gauss):   %.6f Joules\n', trabajo_2_gauss_opt);
fprintf('Trabajo Motor 3 (Simpson): %.6f Joules\n', trabajo_3_simpson_opt);
fprintf('Trabajo Motor 3 (Gauss):   %.6f Joules\n', trabajo_3_gauss_opt);
fprintf('---------------------------------------------------------------------\n\n');

% Gráficas del análisis cinemático para la solución optimizada
figure('Name', 'Análisis Cinemático (Solución Optimizada)');
subplot(2, 2, 1);
plot(tiempo_fino_opt, velocidad_x_opt, 'b-', 'DisplayName', 'Velocidad X');
hold on;
plot(tiempo_fino_opt, velocidad_y_opt, 'r-', 'DisplayName', 'Velocidad Y');
title('Velocidad Cartesiana (Optimizada)');
xlabel('Tiempo (s)'); ylabel('Velocidad (cm/s)'); grid on; legend;

subplot(2, 2, 2);
plot(tiempo_fino_opt, aceleracion_x_opt, 'b-', 'DisplayName', 'Aceleración X');
hold on;
plot(tiempo_fino_opt, aceleracion_y_opt, 'r-', 'DisplayName', 'Aceleración Y');
title('Aceleración Cartesiana (Optimizada)');
xlabel('Tiempo (s)'); ylabel('Aceleración (cm/s^2)'); grid on; legend;

subplot(2, 2, 3);
plot(tiempo_fino_opt, rad2deg(vel_angular_1_opt), 'DisplayName', '\omega_1');
hold on;
plot(tiempo_fino_opt, rad2deg(vel_angular_2_opt), 'DisplayName', '\omega_2');
plot(tiempo_fino_opt, rad2deg(vel_angular_3_opt), 'DisplayName', '\omega_3');
title('Velocidad Angular (Optimizada)');
xlabel('Tiempo (s)'); ylabel('Velocidad Angular (deg/s)'); grid on; legend;

subplot(2, 2, 4);
plot(tiempo_fino_opt, rad2deg(acel_angular_1_opt), 'DisplayName', '\alpha_1');
hold on;
plot(tiempo_fino_opt, rad2deg(acel_angular_2_opt), 'DisplayName', '\alpha_2');
plot(tiempo_fino_opt, rad2deg(acel_angular_3_opt), 'DisplayName', '\alpha_3');
title('Aceleración Angular (Optimizada)');
xlabel('Tiempo (s)'); ylabel('Aceleración Angular (deg/s^2)'); grid on; legend;


% Gráfica del impulso acumulado para la solución optimizada
impulso_acum_x_opt = cumtrapz(tiempo_fino_opt, fuerza_x_opt);
impulso_acum_y_opt = cumtrapz(tiempo_fino_opt, fuerza_y_opt);
magnitud_impulso_acum_opt = sqrt(impulso_acum_x_opt.^2 + impulso_acum_y_opt.^2);

figure('Name', 'Análisis de Impulso (Solución Optimizada)');
plot(tiempo_fino_opt, magnitud_impulso_acum_opt, 'm-', 'LineWidth', 2);
title('Magnitud del Impulso Acumulado vs. Tiempo (Optimizada)');
xlabel('Tiempo (s)');
ylabel('Magnitud del Impulso |J| (N·s)');
grid on;

fprintf('✅ Análisis repetido y gráficos para la solución optimizada completados.\n');