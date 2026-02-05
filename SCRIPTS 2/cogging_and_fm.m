%programma che valuta l'andamento di Fm e cogging force sul mover al
%variare della configurazione. 

%DIFFERENZA TRA MODELLO E DATI SPERIMENTALI
lenght_mover = 150;
pitch_coils = 1;
dim_magnet = [10 15];
diam_coil = 25;
num_coils = 20;
lenght_tot = num_coils*diam_coil + (num_coils-1)*pitch_coils;
pitch = diam_coil + pitch_coils;
dx = 0.1;
m = 0.7; %mass of the cart
F_att = fit(d, F_att1, 'cubicspline');
F_rep = fit(d, F_rep1, 'cubicspline'); 
F_x_dec = zeros(2, length(x_values)); 
F_x_stop = zeros(2, length(x_values)); 
num_magn = 4;
i = 1;
% Create a grid of x and V
x_values = linspace(0, 23, 101); 
V_values = linspace(-24, 24, 100); 
% Define parameters
x_values_plot = 0 : dx : 26;
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
f1_values(x < 13) = F_filtered(x < 13, 1);  % Colonna 1 per x < 13
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

fid = fopen('config salvate.txt', 'r');
if fid == -1
    error('Impossibile aprire il file configurate.txt');
end

best_config_metric = '';
min_metric = Inf;
F_cogging_values = [];
configs = {};

all_var_F_cogging = [];
all_var_tot_fm = [];
all_peak_to_peak_F_cogging = [];

while ~feof(fid)
    % Leggi la riga del file
    line = fgetl(fid);
    
    % Controlla se la riga è vuota
    if ~ischar(line) || isempty(line)
        continue;
    end
    
    % Estrai i valori usando espressioni regolari
    pattern = 'magnete 2: (\d+) // magnete 3: (\d+) // magnete 4: (\d+) // offset: (\d+)';
    tokens = regexp(line, pattern, 'tokens');
    
    if isempty(tokens)
        warning('Formato del file non valido per la riga: %s', line);
        continue;
    end
    
    % Converte i valori in numeri
    values = str2double(tokens{1});
    k = values(1);
    p = values(2);
    b = values(3);
    v = values(4);
    
    % Calcola F_cogging
    F_cogging = compute_force_cogging(k, p, b, v, spline_x, a_opt, b_opt, c_opt, z_opt);
    
    % Calcola varianza e ampiezza picco-picco
    var_F_cogging = var(F_cogging);
    peak_to_peak_F_cogging = max(F_cogging) - min(F_cogging);
    all_var_F_cogging = [all_var_F_cogging, var_F_cogging];
    all_peak_to_peak_F_cogging = [all_peak_to_peak_F_cogging, peak_to_peak_F_cogging];
    
    var_tot_fm = 0;
    for V = 1 : 1 : 24
        Fm_values = compute_Fm(k, p, b, v, alpha_opt, beta_opt, L_opt, k_opt, V);
        var_Fm = var(Fm_values);
        var_tot_fm = var_Fm + var_tot_fm;
    end
    all_var_tot_fm = [all_var_tot_fm, var_tot_fm];
end

fclose(fid);

% Normalizza e calcola la metrica
max_var_fm = max(all_var_tot_fm);
max_var_cogging = max(all_var_F_cogging);
max_peak_to_peak = max(all_peak_to_peak_F_cogging);
fid = fopen('config salvate.txt', 'r');
index = 1;

while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line) || isempty(line)
        continue;
    end
    
    metric = (all_var_tot_fm(index) / max_var_fm) * 0 + (all_var_F_cogging(index) / max_var_cogging) + (all_peak_to_peak_F_cogging(index) / max_peak_to_peak);
    
    if metric < min_metric
        min_metric = metric;
        best_config_metric = line;
        k_save = str2double(regexp(line, '(?<=magnete 2: )\d+', 'match', 'once'));
        p_save = str2double(regexp(line, '(?<=magnete 3: )\d+', 'match', 'once'));
        b_save = str2double(regexp(line, '(?<=magnete 4: )\d+', 'match', 'once'));
        v_save = str2double(regexp(line, '(?<=offset: )\d+', 'match', 'once'));
    end
    index = index + 1;
end

fclose(fid);

% Salva la configurazione migliore in un file
fid_out = fopen('Cogging_e_Fm_ottimizzati.txt', 'w');
if fid_out == -1
    error('Impossibile creare il file best_configurations.txt');
end
fprintf(fid_out, 'Configurazione ottimale secondo la metrica:\n%s\n', best_config_metric);
fclose(fid_out);

% Plot della configurazione migliore
figure;
hold on;
color_idx = 1;
for V = 0 : 4 : 24
    Fm_tot = compute_Fm(k_save, p_save, b_save, v_save, alpha_opt, beta_opt, L_opt, k_opt, V);
    plot(x_values_plot, Fm_tot, 'LineWidth', 2);
    color_idx = color_idx + 1;
end
grid on;
box on;
xlabel('X-position');
ylabel('Controllable force for different voltages');
title(['Best Config - Pitch ', num2str(k_save), ', ', num2str(p_save), ', ', num2str(b_save), ', Offset ', num2str(v_save)]);
hold off;

F_cogging_Fm = compute_force_cogging(k_save, p_save, b_save, v_save, spline_x, a_opt, b_opt, c_opt, z_opt);

figure;
plot(x_values_plot, F_cogging_Fm, 'b', 'LineWidth', 2); % Blu per la configurazione ottimale
grid on;
box on;
xlabel('X-position');
ylabel('Force applied to the system');
title('Best Config - Force Cogging');




function F_tot_cogging = compute_force_cogging(k, p, b, v, spline_x, a_opt, b_opt, c_opt, z_opt)
    pitch_coils = 1;
    diam_coil = 25;
    num_coils = 20;
    pitch = diam_coil + pitch_coils;
    
    dx = 0.1;
    x_values = 0 : dx : 26;
    F_x_acc_modello = zeros(2, length(x_values));
    mover = zeros(2, 4);
    mover(1, :) = [0, k, p, b];
    mover(2, :) = [v, k + v, p + v, b + v];
    
    coils = pitch * (0:num_coils);
    o = 1;
    for x = x_values
         for i = 1:2
            distance = zeros(4, num_coils);
            for s = 1:num_coils  
                for j = 1:4
                    distance(j, s) = coils(s) - mover(i, j);
                end
            end
            
            total_force_modello = 0;
            for j = 1:4
                for s = 1:num_coils
                    if abs(distance(j, s)) < 23
                            if distance(j, s) < 0
                                Fx_modello = -f_cogging(a_opt, b_opt, c_opt, abs(distance(j,s)), 0, z_opt, spline_x);
                            else
                                Fx_modello = f_cogging(a_opt, b_opt, c_opt, abs(distance(j,s)), 0, z_opt, spline_x);
                            end
                        total_force_modello = total_force_modello + Fx_modello;
                    end
                end
            end          
            F_x_acc_modello(i, o) = total_force_modello;
        end
        o = o + 1;
        mover = mover + dx;
    end
    F_tot_cogging = F_x_acc_modello(1, :) + F_x_acc_modello(2, :);
end

function Fm_tot = compute_Fm(k, p, b, v, alpha_opt, beta_opt, L_opt, k_opt, V)
    pitch_coils = 1;
    diam_coil = 25;
    num_coils = 20;
    pitch = diam_coil + pitch_coils;
    dx = 0.1;
    x_values = 0 : dx : 26;
    F_x_acc_modello = zeros(2, length(x_values));
    mover = zeros(2, 4);
    mover(1, :) = [0, k, p, b];
    mover(2, :) = [v, k + v, p + v, b + v]; 
    coils = pitch * (0:num_coils);
    o = 1;
    mode = zeros(1, num_coils);

    for x = x_values
         for i = 1:2
            distance = zeros(4, num_coils);
            for s = 1:num_coils  
                for j = 1:4
                    distance(j, s) = coils(s) - mover(i, j);
                end
            end

            for s = 1:num_coils
                F_total_rep = 0;
                F_total_att = 0;
                for j = 1:4
                    if abs(distance(j, s)) < 23
                        if distance(j, s) < 0
                            F_total_rep = F_total_rep - Fm(alpha_opt, beta_opt, abs(distance(j,s)), -V, L_opt, k_opt);
                            F_total_att = F_total_att - Fm(alpha_opt, beta_opt, abs(distance(j,s)), V, L_opt, k_opt);
                        else
                            F_total_rep = F_total_rep + Fm(alpha_opt, beta_opt, abs(distance(j,s)), -V, L_opt, k_opt);
                            F_total_att = F_total_att + Fm(alpha_opt, beta_opt, abs(distance(j,s)), V, L_opt, k_opt);
                        end
                    end
                end
                if F_total_att > F_total_rep
                    mode(s) = 1;
                else
                    mode(s) = -1;
                end
            end
                
            total_force_modello = 0;

            for j = 1:4
                for s = 1:num_coils
                    if abs(distance(j, s)) < 23
                        if mode(s) < 0
                            if distance(j, s) < 0
                                Fx_modello = -Fm(alpha_opt, beta_opt, abs(distance(j,s)), -V, L_opt, k_opt);
                            else
                                Fx_modello = Fm(alpha_opt, beta_opt, abs(distance(j,s)), -V, L_opt, k_opt);
                            end
                        else
                            if distance(j, s) < 0
                                Fx_modello = -Fm(alpha_opt, beta_opt, abs(distance(j,s)), V, L_opt, k_opt);
                            else
                                Fx_modello = Fm(alpha_opt, beta_opt, abs(distance(j,s)), V, L_opt, k_opt);
                            end
                        end          
                        total_force_modello = total_force_modello + Fx_modello;
                    end
                end
            end
            F_x_acc_modello(i, o) = total_force_modello;
        end
        o = o + 1;
        mover = mover + dx;
    end
    Fm_tot = F_x_acc_modello(1, :) + F_x_acc_modello(2, :);
    
end

function F_cogging = f_cogging(a, b, c, x, V, z, spline_x)
    % Ensure x is a column vector (1xm) and V is a row vector (1xq)
    x = x(:);  % Convert x to column vector
    V = V(:)'; % Convert V to row vector
    
    % Compute f1(x) using the spline
    f1_values = feval(spline_x, x);
    f1_values = f1_values * z;
        
    % Compute f3 over the grid (only for x > 14)
    f3_values = f3(a, b, c, x);
    
    % Expand f1 and f3 to match the dimensions of V
    f1_grid = f1_values * ones(1, length(V));
    f3_grid = f3_values * ones(1, length(V));
    
    % Compute the final function F(x, V)
    F_cogging = (x <= 14) .* (f1_grid) + (x > 14) .* f3_grid;
end

function Fm_val = Fm(alpha, beta, x, V, L, k)
    % Compute f2(x, V) over a grid of x (column) and V (row)
    Fm_val = bsxfun(@times, V, (alpha * exp(-x * beta) + k) .* (sin((pi/2) * (x / L))).^2 );
    Fm_val = (x <= 14) .* Fm_val;
end

function f3_val = f3(a, b, c, x)
    % Compute f3(x) (independent of V)
    f3_val = a * exp(-(((x - b) / c).^2));
end
