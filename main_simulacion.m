%% Simulación del Brazo Robótico (3R)

% 1. Cargar trayectoria
try
    pts = readmatrix('Pichu_trayectoria_all.csv');
catch
    error('Falta el CSV. Ejecuta trazos_por_parte.m primero.');
end

% Escalar al espacio del robot (13 cm)
reach = 13; 
valid = pts(~isnan(pts(:,1)), :);
xmin = min(valid(:,1)); xmax = max(valid(:,1));
ymin = min(valid(:,2)); ymax = max(valid(:,2));

scale = reach / max(xmax-xmin, ymax-ymin);
x = (pts(:,1) - (xmin+xmax)/2) * scale;
y = -(pts(:,2) - (ymin+ymax)/2) * scale;
traj = [x, y];

% 2. Resolver Cinemática Inversa
L1 = 6.0; L2 = 5.5; L3 = 3.0;
th = deg2rad([10; 0; 0]);
hist = cell(size(traj, 1), 1);

for i = 1:size(traj, 1)
    p = traj(i, :);
    if isnan(p(1))
        hist{i} = NaN(3,1);
        continue;
    end
    th = inverse_kinematics_solver(p', th, L1, L2, L3);
    hist{i} = th;
end

% 3. Animación
figure('Name', 'Simulación Brazo');
hold on; axis equal; grid on;
h_brazo = plot(NaN, NaN, 'bo-', 'LineWidth', 2);
h_path = plot(NaN, NaN, 'k-', 'LineWidth', 2);

for i = 1:numel(hist)
    th_i = hist{i};
    if all(isnan(th_i))
        set(h_path, 'XData', [get(h_path, 'XData'), NaN], 'YData', [get(h_path, 'YData'), NaN]);
        continue;
    end
    [x0,y0,x1,y1,x2,y2,x3,y3] = calcular_pos_brazo(th_i, L1, L2, L3);
    set(h_brazo, 'XData', [x0, x1, x2, x3], 'YData', [y0, y1, y2, y3]);
    set(h_path, 'XData', [get(h_path, 'XData'), x3], 'YData', [get(h_path, 'YData'), y3]);
    drawnow; pause(0.01);
end
