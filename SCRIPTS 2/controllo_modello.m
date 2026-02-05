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


%calcolo della forze 
for k = 39 : 1 : 39  
    for p = 72 : 1 : 72
        for b = 110 : 1 : 110
            for v = 10 : 1 : 10
            [flag, F_tot_acc, F_acc_modello] = compute_force_delta_acc(k, p, b, v, F_att, F_rep, spline_x, alpha_opt, beta_opt, a_opt, b_opt, c_opt, L_opt, k_opt, z_opt);
            if flag == 1
                continue;
            end
            [flag, F_tot_dec, F_dec_modello] = compute_force_delta_dec(k, p, b, v, F_att, F_rep, spline_x, alpha_opt, beta_opt, a_opt, b_opt, c_opt, L_opt, k_opt, z_opt);
            if flag == 1
                continue;
            end
            [flag, F_tot_stop, F_stop_modello] = compute_force_delta_stop(k, p, b, v, F_att, F_rep, spline_x, alpha_opt, beta_opt, a_opt, b_opt, c_opt, L_opt, k_opt, z_opt);
            if flag == 1
                continue;
            end            

            figure;
            plot(x_values_plot, F_tot_dec, 'k', 'LineWidth', 2);
            grid on;
            box on;
            xlabel('X-position');
            ylabel('Force applied to the system');
            title(['Force DECELERATION with pitch ', num2str(k), ' and ', num2str(p), ' and ', num2str(b), ' and offset ', num2str(v)]);
            hold on;
            plot(x_values_plot, F_dec_modello, 'g', 'LineWidth', 2);


            figure;
            plot(x_values_plot, F_tot_stop, 'k', 'LineWidth', 2);
            grid on;
            box on;
            xlabel('X-position');
            ylabel('Force applied to the system');
            title(['Force STOP with pitch ', num2str(k), ' and ', num2str(p), ' and ', num2str(b), ' and offset ', num2str(v)]);
            hold on;
            plot(x_values_plot, F_stop_modello, 'g', 'LineWidth', 2);

            figure;
            plot(x_values_plot, F_tot_acc, 'k', 'LineWidth', 2);
            grid on;
            box on;
            xlabel('X-position');
            ylabel('Force applied to the system');
            title(['Force ACCELERATION with pitch ', num2str(k), ' and ', num2str(p), ' and ', num2str(b), ' and offset ', num2str(v)]);
            hold on;
            plot(x_values_plot, F_acc_modello, 'g', 'LineWidth', 2);
            

            prova = (F_tot_acc - F_cogging_Fm);

            figure;
            plot(x_values_plot, prova, 'k', 'LineWidth', 2);
            grid on;
            box on;
            xlabel('X-position');
            ylabel('Force applied to the system');
            title('Expected controllable force(green) vs real one');
            hold on;
            plot(x_values_plot, Fm_tot, 'g', 'LineWidth', 2);






            fprintf('configurazione da salvare: magnete 2: %d // magnete 3: %d // magnete 4: %d // offset: %d \n', k, p, b, v);

            end
        end
    end
end




function [flag, F_tot_acc, F_tot_acc_modello] = compute_force_delta_acc(k, p, b, v, F_att, F_rep, spline_x, alpha_opt, beta_opt, a_opt, b_opt, c_opt, L_opt, k_opt, z_opt)
    pitch_coils = 1;
    diam_coil = 25;
    num_coils = 20;
    pitch = diam_coil + pitch_coils;
    
    dx = 0.1;
    x_values = 0 : dx : 26;
    F_x_acc = zeros(2, length(x_values));
    F_x_acc_modello = zeros(2, length(x_values));
    F_tot_acc = zeros(1, length(x_values));
    F_acc_modello = zeros(1, length(x_values));
    mover = zeros(2, 4);
    mover(1, :) = [0, k, p, b];
    mover(2, :) = [v, k + v, p + v, b + v];
    
    coils = pitch * (0:num_coils);
    o = 1;
    flag = 0;
    for x = x_values
         for i = 1:2
            distance = zeros(4, num_coils);
            mode = zeros(1, num_coils);
            
            for s = 1:num_coils  
                for j = 1:4
                    distance(j, s) = coils(s) - mover(i, j);
                end
            end
            
            total_force = 0;
            total_force_modello = 0;

            for s = 1:num_coils
                F_total_rep = 0;
                F_total_att = 0;
                
                for j = 1:4
                    if abs(distance(j, s)) < 23
                        if distance(j, s) < 0
                            F_total_rep = F_total_rep - F_rep(abs(distance(j, s)));
                            F_total_att = F_total_att - F_att(abs(distance(j, s)));
                        else
                            F_total_rep = F_total_rep + F_rep(abs(distance(j, s)));
                            F_total_att = F_total_att + F_att(abs(distance(j, s)));
                        end
                    end
                end
                
                if F_total_att > F_total_rep
                    mode(s) = 1;
                else
                    mode(s) = -1;
                end
            end
            
            for j = 1:4
                for s = 1:num_coils
                    if abs(distance(j, s)) < 23
                        if mode(s) < 0
                            if distance(j, s) < 0
                                Fx = -F_rep(abs(distance(j, s)));
                                Fx_modello = -f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, abs(distance(j,s)), -24, L_opt, k_opt, z_opt, spline_x);
                            else
                                Fx = F_rep(abs(distance(j, s)));
                                Fx_modello = f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, abs(distance(j,s)), -24, L_opt, k_opt, z_opt, spline_x);
                            end
                        else
                            if distance(j, s) < 0
                                Fx = -F_att(abs(distance(j, s)));
                                Fx_modello = -f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, abs(distance(j,s)), 24, L_opt, k_opt, z_opt, spline_x);
                            else
                                Fx = F_att(abs(distance(j, s)));
                                Fx_modello = f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, abs(distance(j,s)), 24, L_opt, k_opt, z_opt, spline_x);
                            end
                        end
                        total_force = total_force + Fx;
                        total_force_modello = total_force_modello + Fx_modello;
                    end
                end
            end          
            F_x_acc(i, o) = total_force;
            F_x_acc_modello(i, o) = total_force_modello;
        end
        o = o + 1;
        mover = mover + dx;
    end
    F_tot_acc = F_x_acc(1, :) + F_x_acc(2, :);
    F_tot_acc_modello = F_x_acc_modello(1, :) + F_x_acc_modello(2, :);
end

function [flag, F_tot_dec, F_tot_dec_modello] = compute_force_delta_dec(k, p, b, v, F_att, F_rep, spline_x, alpha_opt, beta_opt, a_opt, b_opt, c_opt, L_opt, k_opt, z_opt)
    pitch_coils = 1;
    diam_coil = 25;
    num_coils = 20;
    pitch = diam_coil + pitch_coils;
    dx = 0.1;
    x_values = 0:dx:26;
    F_x_acc_modello = zeros(2, length(x_values));
    F_x_acc = zeros(2, length(x_values));
    F_tot_dec_modello = zeros(1, length(x_values));
    F_tot_dec = zeros(1, length(x_values));

    mover = zeros(2, 4);
    mover(1, :) = [0, k, p, b];
    mover(2, :) = [v, k + v, p + v, b + v];

    coils = pitch * (0:num_coils);
    o = 1;
    flag = 0;

     for x = x_values
         for i = 1:2
            distance = zeros(4, num_coils);
            mode = zeros(1, num_coils);

            for s = 1:num_coils  
                for j = 1:4
                    distance(j, s) = coils(s) - mover(i, j);
                end
            end

            total_force = 0;
            total_force_modello = 0;

            for s = 1:num_coils
                F_total_rep = 0;
                F_total_att = 0;

                for j = 1:4
                    if abs(distance(j, s)) < 23
                        if distance(j, s) < 0
                            F_total_rep = F_total_rep - F_rep(abs(distance(j, s)));
                            F_total_att = F_total_att - F_att(abs(distance(j, s)));
                        else
                            F_total_rep = F_total_rep + F_rep(abs(distance(j, s)));
                            F_total_att = F_total_att + F_att(abs(distance(j, s)));
                        end
                    end
                end

                if F_total_att > F_total_rep
                    mode(s) = -1;
                else
                    mode(s) = 1;
                end
            end

            for j = 1:4
                for s = 1:num_coils
                    if abs(distance(j, s)) < 23
                        if mode(s) < 0
                            if distance(j, s) < 0
                                Fx = -F_rep(abs(distance(j, s)));
                                Fx_modello = -f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, abs(distance(j,s)), -24, L_opt, k_opt, z_opt, spline_x);
                            else
                                Fx = F_rep(abs(distance(j, s)));
                                Fx_modello = f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, abs(distance(j,s)), -24, L_opt, k_opt, z_opt, spline_x);
                            end
                        else
                            if distance(j, s) < 0
                                Fx = -F_att(abs(distance(j, s)));
                                Fx_modello = -f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, abs(distance(j,s)), 24, L_opt, k_opt, z_opt, spline_x);
                            else
                                Fx = F_att(abs(distance(j, s)));
                                Fx_modello = f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, abs(distance(j,s)), 24, L_opt, k_opt, z_opt, spline_x);
                            end
                        end
                        total_force = total_force + Fx;
                        total_force_modello = total_force_modello + Fx_modello;
                    end
                end
            end          
            F_x_acc(i, o) = total_force;
            F_x_acc_modello(i, o) = total_force_modello;
        end
        o = o + 1;
        mover = mover + dx;
    end
    F_tot_dec = F_x_acc(1, :) + F_x_acc(2, :);
    F_tot_dec_modello = F_x_acc_modello(1, :) + F_x_acc_modello(2, :);
end

function [flag, F_tot_stop, F_tot_stop_modello] = compute_force_delta_stop(k, p, b, v, F_att, F_rep, spline_x, alpha_opt, beta_opt, a_opt, b_opt, c_opt, L_opt, k_opt, z_opt)
    pitch_coils = 1;
    diam_coil = 25;
    num_coils = 20;
    pitch = diam_coil + pitch_coils;
    dx = 0.1;
    x_values = 0:dx:26;
    F_x_acc = zeros(2, length(x_values));
    F_tot_stop = zeros(1, length(x_values));
    F_tot_stop_modello = zeros(1, length(x_values));

    mover = zeros(2, 4);
    mover(1, :) = [0, k, p, b];
    mover(2, :) = [v, k + v, p + v, b + v];

    coils = pitch * (0:num_coils);
    flag = 0;
    o = 1;
    for x = x_values
         for i = 1:2
            distance = zeros(4, num_coils);
            mode = zeros(1, num_coils);

            for s = 1:num_coils  
                for j = 1:4
                    distance(j, s) = coils(s) - mover(i, j);
                end
            end

            total_force = 0;
            total_force_modello = 0;

            for s = 1:num_coils
                F_total_rep = 0;
                F_total_att = 0;

                for j = 1:4
                    if abs(distance(j, s)) < 23
                        if distance(j, s) < 0
                            F_total_rep = F_total_rep - F_rep(abs(distance(j, s)));
                            F_total_att = F_total_att - F_att(abs(distance(j, s)));
                        else
                            F_total_rep = F_total_rep + F_rep(abs(distance(j, s)));
                            F_total_att = F_total_att + F_att(abs(distance(j, s)));
                        end
                    end
                end

                if abs(F_total_att) > abs(F_total_rep)
                    mode(s) = -1;
                else
                    mode(s) = 1;
                end
            end
            for j = 1:4
                for s = 1:num_coils
                    if abs(distance(j, s)) < 23
                        if mode(s) < 0
                            if distance(j, s) < 0
                                Fx = -F_rep(abs(distance(j, s)));
                                Fx_modello = -f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, abs(distance(j,s)), 0, L_opt, k_opt, z_opt, spline_x);
                            else
                                Fx = F_rep(abs(distance(j, s)));
                                Fx_modello = f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, abs(distance(j,s)), 0, L_opt, k_opt, z_opt, spline_x);
                            end
                        else
                            if distance(j, s) < 0
                                Fx = -F_att(abs(distance(j, s)));
                                Fx_modello = -f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, abs(distance(j,s)), 0, L_opt, k_opt, z_opt, spline_x);
                            else
                                Fx = F_att(abs(distance(j, s)));
                                Fx_modello = f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, abs(distance(j,s)), 0, L_opt, k_opt, z_opt, spline_x);
                            end
                        end
                        total_force = total_force + Fx;
                        total_force_modello = total_force_modello + Fx_modello;
                    end
                end
            end          
            F_x_acc(i, o) = total_force;
            F_x_acc_modello(i, o) = total_force_modello;
        end
        o = o + 1;
        mover = mover + dx;
    end
    F_tot_stop = F_x_acc(1, :) + F_x_acc(2, :);
    F_tot_stop_modello = F_x_acc_modello(1, :) + F_x_acc_modello(2, :);
end       

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
