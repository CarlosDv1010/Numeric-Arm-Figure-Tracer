function trazos_por_parte()
%% ===== Parámetros =====
imgPath       = 'fig5.jpg';  % <- tu ruta
minPts        = 40;          % descarta trazos muy cortos
doResample    = true;       % true para remuestrear cada trazo a Nfijo
NperTrace     = 300;         % usado solo si doResample=true
cannyThresh   = [0.015 0.28];% [low high] para edge('canny'); ajusta si hace falta
gaussSigma    = 0.6;         % suavizado suave antes de Canny
csvPrefix     = 'Pichu';     % prefijo para los CSV resultantes

%% ===== Cargar imagen =====
assert(exist(imgPath,'file')==2,'No se encontró: %s',imgPath);
try
    [RGB,~,A] = imread(imgPath);
catch
    RGB = imread(imgPath); A = [];
end
RGB = im2double(RGB);
if size(RGB,3)==1, RGB = repmat(RGB,[1 1 3]); end
Igray = rgb2gray(RGB);
[H,W,~] = size(RGB); %#ok<ASGLU>

figure('Name','Imagen original'); imshow(RGB); title('Imagen original');

%% ===== Bordes Canny (LA figura "buena") =====
Isoft  = imgaussfilt(Igray, gaussSigma);
BWedge = edge(Isoft, 'canny', cannyThresh);

% Si hay alpha o fondo muy distinto, recorta por máscara (no altera bordes)
if ~isempty(A)
    BWedge = BWedge & (A>0);
else
    % recorte suave por color de esquinas (opcional; comenta si no hace falta)
    corners = [RGB(1,1,:); RGB(1,end,:); RGB(end,1,:); RGB(end,end,:)];
    bg = squeeze(median(corners,1))';
    dist = sqrt(sum((reshape(RGB,[],3)-bg).^2,2));
    fgMask = reshape(dist > 0.05, size(RGB,1), size(RGB,2));
    BWedge = BWedge & fgMask;
end

figure('Name','Bordes (usar estos para trayectorias)'); imshow(BWedge);
title('Bordes Canny: estas curvas = trayectorias');

%% ===== Trajos = Contornos de los bordes (tal cual) =====
B = bwboundaries(BWedge, 8, 'noholes'); % cell de contornos en [fila col]
if isempty(B), error('No se detectaron bordes. Ajusta cannyThresh/gaussSigma.'); end

% Filtra por tamaño mínimo
lens = cellfun(@(c) size(c,1), B);
keep = lens >= minPts;
B    = B(keep);
lens = lens(keep);

% Orden de trazos (opcional): largo a corto
[~,ord] = sort(lens,'descend'); B = B(ord); lens = lens(ord);

% Vista rápida
figure('Name','Trazos extraídos de los bordes'); imshow(BWedge); hold on;
for j=1:numel(B), C=B{j}; plot(C(:,2),C(:,1),'.'); end
title(sprintf('Trazos: %d', numel(B))); hold off;

%% ===== Armar y (opcional) remuestrear cada trazo =====
trayectorias_img = cell(numel(B),1);
for j=1:numel(B)
    C = B{j};                    % [fila(y), col(x)]
    xy = [double(C(:,2)), double(C(:,1))];

    % quitar duplicados consecutivos mínimos
    dup = [false; hypot(diff(xy(:,1)), diff(xy(:,2))) < 1e-12];
    xy  = xy(~dup,:);

    if doResample
        closed = hypot(xy(1,1)-xy(end,1), xy(1,2)-xy(end,2)) < 1.0;
        xy = resample_polyline(xy, NperTrace, closed);
    end
    trayectorias_img{j} = xy;   % *** puntos EXACTOS de los bordes ***
end

%% ===== Guardado: un CSV por trazo + uno único con NaN =====
for j=1:numel(trayectorias_img)
    fname = sprintf('%s_T%d.csv', csvPrefix, j);
    writematrix(trayectorias_img{j}, fname);  % columnas: x, y (pixeles)
end

trayectoria_all = [];
for j=1:numel(trayectorias_img)
    trayectoria_all = [trayectoria_all; trayectorias_img{j}; [NaN NaN]]; %#ok<AGROW>
end
writematrix(trayectoria_all, sprintf('%s_trayectoria_all.csv', csvPrefix));

fprintf('\n✅ Guardado CSVs en coordenadas de imagen (pixeles):\n');
fprintf('  - %s_T#.csv (un archivo por trazo)\n', csvPrefix);
fprintf('  - %s_trayectoria_all.csv (con NaN separadores)\n\n', csvPrefix);

%% ===== (Opcional) Convertir a marco del robot =====
% Si luego quieres el marco del robot (escala/centro), descomenta:
%{
reach_target = 12;  % largo máximo deseado en unidades del robot
allP = vertcat(trayectorias_img{:});
xmin=min(allP(:,1)); xmax=max(allP(:,1));
ymin=min(allP(:,2)); ymax=max(allP(:,2));
w=xmax-xmin; h=ymax-ymin; xmid=(xmax+xmin)/2; ymid=(ymax+ymin)/2;
scale = reach_target / max(w,h);

trayectorias_robot = cellfun(@(P) scale*[P(:,1)-xmid, -(P(:,2)-ymid)], ...
                             trayectorias_img, 'UniformOutput', false);
save('trayectorias_robot.mat','trayectorias_robot');
%}
end

%% ===== Helper =====
function xy2 = resample_polyline(xy, N, closed)
    if nargin<3, closed=false; end
    if closed, xy=[xy; xy(1,:)]; end
    seg = sqrt(sum(diff(xy,1,1).^2,2));
    s   = [0; cumsum(seg)];
    if s(end)==0, xy2=repmat(xy(1,:),N,1); return; end
    
    % --- INICIO DE LA CORRECCIÓN ---
    % Identificar y eliminar puntos donde la distancia acumulada 's' no aumenta.
    % Esto elimina los puntos duplicados que causan el error en interp1.
    [s_unique, unique_indices, ~] = unique(s, 'stable');
    
    % Si todos los puntos eran duplicados excepto el primero, manejamos el caso.
    if numel(s_unique) < 2
        xy2 = repmat(xy(1,:), N, 1);
        return;
    end
    
    % Nos quedamos solo con los puntos únicos
    xy_unique = xy(unique_indices, :);
    s = s_unique; % 's' ahora es estrictamente creciente
    xy = xy_unique;
    % --- FIN DE LA CORRECCIÓN ---

    sN = linspace(0, s(end), N).';
    xN = interp1(s, xy(:,1), sN, 'linear');
    yN = interp1(s, xy(:,2), sN, 'linear');
    xy2 = [xN, yN];
end