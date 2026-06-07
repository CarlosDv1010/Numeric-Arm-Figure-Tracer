function trazos_por_parte()
% Extrae trayectorias de una imagen usando Canny

imgPath = 'fig5.jpg'; 
minPts = 40;          
NperTrace = 300;         
csvPrefix = 'Pichu';     

img = im2double(imread(imgPath));
if size(img,3)==1, img = repmat(img,[1 1 3]); end
Igray = rgb2gray(img);

% Detectar bordes
Isoft = imgaussfilt(Igray, 0.6);
BW = edge(Isoft, 'canny', [0.015 0.28]);

% Obtener contornos
B = bwboundaries(BW, 8, 'noholes');
lens = cellfun(@(c) size(c,1), B);
B = B(lens >= minPts);

% Procesar y remuestrear
traj = cell(numel(B),1);
for j=1:numel(B)
    C = B{j};
    xy = [double(C(:,2)), double(C(:,1))];
    % Quitar duplicados
    dup = [false; hypot(diff(xy(:,1)), diff(xy(:,2))) < 1e-12];
    xy = xy(~dup,:);
    
    % Remuestrear para suavidad
    closed = hypot(xy(1,1)-xy(end,1), xy(1,2)-xy(end,2)) < 1.0;
    traj{j} = resample_polyline(xy, NperTrace, closed);
end

% Guardar CSV único con separadores NaN
all_traj = [];
for j=1:numel(traj)
    all_traj = [all_traj; traj{j}; [NaN NaN]]; 
end
writematrix(all_traj, sprintf('%s_trayectoria_all.csv', csvPrefix));

end

function xy2 = resample_polyline(xy, N, closed)
    if closed, xy=[xy; xy(1,:)]; end
    seg = sqrt(sum(diff(xy,1,1).^2,2));
    s = [0; cumsum(seg)];
    if s(end)==0, xy2=repmat(xy(1,:),N,1); return; end
    
    [s_unq, idx_unq] = unique(s, 'stable');
    if numel(s_unq) < 2
        xy2 = repmat(xy(1,:), N, 1);
        return;
    end
    
    sn = linspace(0, s_unq(end), N).';
    xy2 = [interp1(s_unq, xy(idx_unq,1), sn), interp1(s_unq, xy(idx_unq,2), sn)];
end
