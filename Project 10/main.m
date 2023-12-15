clear all;
clc;

% Charger les données à partir du fichier Alicona
addpath("AliconaDataHandling\")
data = AliconaReader()


% Vérifier si les données de texture sont présentes( les données detextures dans le fichier data sont trés grand, j'ai essais de prendre une  sous_matrice pour mettre une petite quantité dedans)
if isfield(data, 'TextureData')
    % Extraire une sous-matrice des données de texture (par exemple, les 100 premières lignes et colonnes)
    subset_texture = data.TextureData(1:1354, 1:1354, :);

    % Afficher la sous-matrice des données de texture avec imshow
    imshow(subset_texture);
    title('Subset of Texture Data');
else
    disp('Les données de texture ne sont pas présentes dans la structure.');
end


%% convertir les données de DepthData en mm et faire sortir les x y
clc;
depthMatrix = data.DepthData;

% Obtenez les coordonnées x et y
[x, y] = meshgrid(1:size(depthMatrix, 2), 1:size(depthMatrix, 1));

pixelSizeX = 0.00000199094008884966; 
pixelSizeY = 0.00000199094008884966;

% Convertir les indices en millimètres en utilisant l'échelle arbitraire
x_mm = x * pixelSizeX;
y_mm = y * pixelSizeY;
% Obtenez les valeurs de profondeur en millimètres

z_mm = reshape(depthMatrix .* 0.00000197222343385677 + 0.00000000765763348575, size(depthMatrix)); %<depthQuality minVal="0.00000000000000000000" filterVal="0.00000197222343385677" shift="0.00000000765763348575"/>
% on travaille tjrs de ligne 6 à 1836 pour éviter les erreur de bord
x = x_mm(6:1800, 6:1800);
y = y_mm(6:1800, 6:1800);
z = z_mm(6:1800, 6:1800);
%%
% Afficher la topographie de surface avec les indicateurs de rugosité
AliconaPlot(data);


%Afficher la légende avec les indicateurs de rugosité
title(sprintf('Rugosité de surface - %s', data.Header.Type));

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');

%% 2. Dégauchissage du nuage de points (utilisation d’un plan des moindres carrés).

% Obtenir la description des données
description = GetAliconaDataDescription(data);
disp(description);

% Visualisation des résultats
figure;
scatter3(x_mm(6:1800,:), y_mm(6:1800,:), z_mm(6:1800,:), 1, 'b.');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Profondeur (mm)');
title('Nuage de points 3D');
grid on;
hold on;

%% Tracer le plan de moindre carrées
clc;
mdl = fitlm([x(:), y(:)], z(:), 'linear');% fct de moindre carrées....
%methode basique:

% Calculer les moyennes
mean_x = mean(x);
mean_y = mean(y);
mean_z = mean(z);

% Calculer les coefficients du plan des moindres carrés
a = sum((x - mean_x) .* (y - mean_y) .* (z - mean_z)) / sum((x - mean_x).^2);
b = sum((x - mean_x) .* (y - mean_y) .* (z - mean_z)) / sum((y - mean_y).^2);
c = mean_z - a * mean_x - b * mean_y;

% Créer une grille pour le plan
[xGrid, yGrid] = meshgrid(linspace(min(x(:)), max(x(:)), 100), linspace(min(y(:)), max(y(:)), 100));

% Calculer les valeurs de z sur la grille en utilisant l'équation du plan
zGrid = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xGrid + mdl.Coefficients.Estimate(3) * yGrid;


% Visualiser le nuage de points et le plan ajusté
figure;
scatter3(x(:), y(:), z(:), 1, 'b.');
hold on;

mesh(xGrid, yGrid, zGrid, 'EdgeColor', 'r', 'FaceAlpha', 0.5);% Tracer le plan ajusté

xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title('Nuage de points et plan des moindres carrés');
legend('Nuage de points', 'Plan des moindres carrés', 'Location', 'Best');
grid on;
hold off;

%%

% Modèle de plan des moindres carrés
mdl = fitlm([x(:), y(:)], z(:), 'linear');

% Prédire les valeurs du modèle (plan des moindres carrés)
z_pred = predict(mdl, [x(:), y(:)]);

% Reshape z_pred pour obtenir la même forme que z
z_pred_reshaped = reshape(z_pred, size(z));

% Dégauchir les données
z_detrended = z - z_pred_reshaped;

% Tracer les résultats
figure;

% Tracer le nuage de points original
subplot(2, 1, 1);
scatter3(x(:), y(:), z(:), 1, 'filled');
title('Données Originales');
xlabel('X');
ylabel('Y');
zlabel('Z');

% Tracer le nuage de points dégauchi
subplot(2, 1, 2);
scatter3(x(:), y(:), z_detrended(:), 1, 'filled');
title('Données Dégauchies');
xlabel('X');
ylabel('Y');
zlabel('Dégauchissement de Z');
%% Critères de rugosité surfacique.
clc;

% Critères de rugosité surfacique.



% Calculer les critères de rugosité surfacique
Sa = (1 / numel(z)) * sum(abs(z - mean(z)), 'all');
Sz = max(z, [], 'all') + abs(min(z, [], 'all'));
Sq = std(z, 0, 'all');  % Écart-type des hauteurs (Moyenne quadratique hauteur)

% Utiliser la Toolbox Statistics pour calculer la curtose et l'asymétrie
if license('test', 'Statistics_Toolbox')
    Ssk = skewness(z(:), 0);  % skewness avec biais = 0
    Sku = kurtosis(z(:), 0);  % kurtosis avec biais = 0
else
    disp('La Toolbox Statistics de MATLAB est requise pour calculer Ssk et Sku.');
    Ssk = NaN;
    Sku = NaN;
end

% Calculer les critères additionnels
Sp = max(z, [], 'all');  % Hauteur de crête maximale
Sv = abs(min(z, [], 'all'));  % Hauteur de puits maximale

% Afficher les résultats
disp('Critères de rugosité surfacique:');
disp(['Sa : ' num2str(Sa)]);
disp(['Sz : ' num2str(Sz)]);
disp(['Sq : ' num2str(Sq)]);
disp(['Ssk : ' num2str(Ssk)]);
disp(['Sku : ' num2str(Sku)]);
disp(['Sp : ' num2str(Sp)]);
disp(['Sv : ' num2str(Sv)]);

 %%  FILTAGE 

R = 150e-6; % Rayon du disque (50 micromètres)
xd = 2.2297e-6; % Pas de discrétisation en X
yd = xd; % Pas de discrétisation en Y (supposé égal à xd)
% Dimensions de la matrice de surface
[Ny, Nx] = size(z_detrended);

% Calculer le nombre de lectures du diamètre du disque en X et Y
Mx = round(2 * R / xd);
My = Mx; % Utiliser la même valeur pour Y

% Préparer une matrice pour stocker les résultats du filtrage
b = zeros(size(z_detrended));

% Assurer que Mx et My sont impairs pour avoir un centre bien défini
Mx = Mx + mod(Mx, 2) - 1;
My = My + mod(My, 2) - 1;

% Boucle de filtrage morphologique
for nx = 1 + floor(Mx/2) : Nx - ceil(Mx/2)
    for ny = 1 + floor(My/2) : Ny - ceil(My/2)
        h = -inf(Mx, My);  % Initialiser h à -inf pour chaque point de la surface

        for mx = 1:Mx
            for my = 1:My
                idx_x = nx - ceil(Mx/2) + mx;
                idx_y = ny - ceil(My/2) + my;

                % Calculer la distance carrée du centre du disque
                distSquared = (xd * (-mx + ceil(Mx/2)))^2 + (yd * (-my + ceil(My/2)))^2;

                if distSquared <= R^2
                    % Calculer la hauteur si à l'intérieur du rayon du disque
                    deltaX = x( idx_x )-x(nx);
                    deltaY = y( idx_y )-y(ny);
                    h(mx, my) = z_detrended(idx_y, idx_x) + sqrt(R^2 - deltaX^2 - deltaY^2) - R;
                end
            end
        end

        % Trouver la hauteur maximale pour les différentes lectures
        b(ny, nx) = max(h(:));
    end
end


% Créer une figure
figure;

% Afficher le nuage de points original en bleu
scatter3(x, y, z_detrended, 1, 'filled', 'b');
hold on;

% Afficher le nuage de points filtré en rouge
scatter3(x, y, b, 1, 'filled', 'r');

% Configurer le titre et les étiquettes des axes
title('Comparaison entre Nuage de points Original et dilater');
xlabel('X');
ylabel('Y');
zlabel('Hauteur');

% Légende pour distinguer les ensembles de données
legend('Original', 'dilater', 'Location', 'Best');

% Réglage des paramètres
grid on;
hold off;
% erosion 
   b_eroded = zeros(size(z_detrended));

for nx = 1 + floor(Mx/2) : Nx - ceil(Mx/2)
    for ny = 1 + floor(My/2) : Ny - ceil(My/2)
        h_eroded = +inf(Mx, My);  % Initialiser h à +inf pour chaque point de la surface car on va prendre le min apres

        for mx = 1:Mx
            for my = 1:My
                idx_x = nx - ceil(Mx/2) + mx;
                idx_y = ny - ceil(My/2) + my;

                % Calculer la distance carrée du centre du disque
                distSquared = (xd * (-mx + ceil(Mx/2)))^2 + (yd * (-my + ceil(My/2)))^2;

                if distSquared <= R^2
                    % Calculer la hauteur si à l'intérieur du rayon du disque
                    deltaX = x( idx_x )-x(nx);
                    deltaY = y( idx_y )-y(ny);
                h_eroded(mx, my) =b(idx_y, idx_x) - sqrt(R^2 - deltaX^2 - deltaY^2) + R;

                
                end
            end
        end
        b_eroded(ny, nx) = min( h_eroded(:));
% non_zero_h = h_eroded(h_eroded > 0);  % Exclure les zéros
%         if ~isempty(non_zero_h)
%             b_eroded(ny, nx) = min( h_eroded(:));
%         end
    end
end

%Créer une figure
figure;

% Afficher le nuage de points original en bleu
scatter3(x, y, z_detrended, 1, 'filled', 'b');
hold on;

% Afficher le nuage de points erosionné en rouge
scatter3(x, y, b_eroded, 1, 'filled', 'r');

title('Comparaison entre Nuage de points Original et érosionné');
xlabel('X');
ylabel('Y');
zlabel('Hauteur');

% Légende pour distinguer les ensembles de données
legend('Original', 'érosionné', 'Location', 'Best');

% Réglage des paramètres
grid on;
hold off;

