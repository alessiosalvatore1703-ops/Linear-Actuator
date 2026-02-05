
function F = f(alpha, beta, a, b, c, x, V, L, k, z, spline_x)
    % Ensure x is a column vector (1xm) and V is a row vector (1xq)
    x = x(:);  % Convert x to column vector
    V = V(:)'; % Convert V to row vector
    
    % Compute f1(x) using the spline
    f1_values = feval(spline_x, x);
    f1_values = f1_values * z;
    
    % Compute f2 over the grid
    f2_values = f2(alpha, beta, x, V, L, k);
    
    % Compute f3 over the grid (only for x > 14)
    f3_values = f3(a, b, c, x);
    
    % Expand f1 and f3 to match the dimensions of V
    f1_grid = f1_values * ones(1, length(V));
    f3_grid = f3_values * ones(1, length(V));
    
    % Compute the final function F(x, V)
    F = (x <= 14) .* (f1_grid + f2_values) + (x > 14) .* f3_grid;
end

function f2_val = f2(alpha, beta, x, V, L, k)
    % Compute f2(x, V) over a grid of x (column) and V (row)
    f2_val = bsxfun(@times, V, (alpha * exp(-x * beta) + k) .* (sin((pi/2) * (x / L))).^2 );
end

function f3_val = f3(a, b, c, x)
    % Compute f3(x) (independent of V)
    f3_val = a * exp(-(((x - b) / c).^2));
end

% Create a grid of x and V
x_values = linspace(0, 23, 100); 
V_values = linspace(-24, 24, 100); 
% Define parameters

% Caricamento dati da Excel
filename = 'cubic_pm.xlsx';
data = readmatrix(filename);
x = data(2:end, 1);  % Posizioni x (escludendo l'intestazione)
V = [0, -12, 12, -24, 24];  % Valori di V (esclusa l'intestazione)
F = data(2:end, 2:end);  % Forze (matrice con le forze per i vari V)


% Selezione dei valori di x tra 13 e 15 e calcolo della media
mask_13_15 = (x >= 13) & (x <= 15);
mean_F_13_15 = mean(F(mask_13_15, :), 2);

% Costruzione del vettore dei valori di f1(x)
f1_values = zeros(size(x));
f1_values(x < 13) = F(x < 13, 1);  % Colonna 1 per x < 13
f1_values(mask_13_15) = mean_F_13_15;  % Media per x tra 13 e 15

% Creazione della spline
spline_x = fit(x, f1_values, 'cubicspline');

% Parametri ottimizzati da file force_optimization.m
alpha_opt = optimized_params(1);
beta_opt = optimized_params(2);
a_opt = optimized_params(3);
b_opt = optimized_params(4);
c_opt = optimized_params(5);
L_opt = optimized_params(6);
k_opt = optimized_params(7);
z_opt = optimized_params(8);


% Compute F
F_values = f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, x_values, V_values, L_opt, k_opt, z_opt, spline_x);

% Plot the result
figure;
surf(V_values, x_values, F_values);
xlabel('V'); ylabel('x'); zlabel('F(x, V)');
title('Surface plot of F(x, V)');
colorbar;



