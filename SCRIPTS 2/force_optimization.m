% Caricamento dati da Excel
filename = 'cubic_pm.xlsx';
data = readmatrix(filename);
x = data(2:end, 1);  % Posizioni x (escludendo l'intestazione)
V = [0, -12, 12, -24, 24];  % Valori di V (esclusa l'intestazione)
F_filtered = data(2:end, 2:end);  % Forze (matrice con le forze per i vari V)

% Filtro di Savitzky-Golay per ridurre il rumore nei dati
window_size = 5;  % Imposta la dimensione della finestra
polynomial_order = 3;  % Ordine del polinomio
%F_filtered = sgolayfilt(F, polynomial_order, window_size);

% Selezione dei valori di x tra 13 e 15 e calcolo della media
mask_13_15 = (x >= 13) & (x <= 15);
mean_F_13_15 = mean(F_filtered(mask_13_15, :), 2);

% Costruzione del vettore dei valori di f1(x)
f1_values = zeros(size(x));
f1_values(x < 13) = F_filtered(x < 13, 1);  % Colonna 1 per x < 13
f1_values(mask_13_15) = mean_F_13_15;  % Media per x tra 13 e 15

% Creazione della spline
spline_x = fit(x, f1_values, 'cubicspline');

% Funzione f2(x, V) dipendente da x e V
f2 = @(alpha, beta, x, V, L, k) V .* (alpha *  exp(-x * beta) + k ) .* (sin(pi/2 * (x/L))).^2 ;

% Funzione f3(x) parabola indipendente da V (solo in [15, 23])
f3 = @(a, b, c, x) a * exp(-(((x - b)/c).^2));

% Funzione completa F(x, V) = f1(x) + f2(x, V) + f3(x)
f = @(alpha, beta, a, b, c, x, V, L, k, z) ...
    (x <= 14) .* (ones(24, 1) * ones(1, 5)) .* (z * feval(spline_x, x)*ones(1,5) + f2(alpha, beta, x, V, L, k)) + ...
    (x > 14) .* (ones(24, 1) * ones(1, 5)) .* (f3(a, b, c, x)*ones(1, 5));

% Funzione di errore (MSE) che combina tutte le componenti
mse = @(params, x, V, F) sum((F - f(params(1), params(2), params(3), params(4), params(5), x, V, params(6), params(7), params(8))).^2, "all") / numel(F);


num_trials = 100; % Numero di tentativi
best_params = [];
best_mse = Inf;

for i = 1:num_trials
    % Generazione di nuovi valori iniziali attorno a quelli originali
    initial_params = [-0.1667, 0.005, 7.13558828629857, 16.9812285795709, 2.74673279604822, -0.0273, 1, 1] ...
                        + 10*randn(1, 8); % Aggiunge rumore casuale
    
    % Ottimizzazione
    [params, mse_value] = fminsearch(@(p) mse(p, x, V, F_filtered), initial_params);
    
    % Aggiornamento del miglior risultato
    if mse_value < best_mse
        best_mse = mse_value;
        best_params = params;
    end
end

optimized_params = best_params;



% Parametri ottimizzati
alpha_opt = optimized_params(1);
beta_opt = optimized_params(2);
a_opt = optimized_params(3);
b_opt = optimized_params(4);
c_opt = optimized_params(5);
L_opt = optimized_params(6);
k_opt = optimized_params(7);
z_opt = optimized_params(8);

% Calcolare la forza predetta con i parametri ottimizzati
% Ora possiamo calcolare F_pred senza cicli, direttamente con matrici
% Calcolare le forze predette
F_pred = f(alpha_opt, beta_opt, a_opt, b_opt, c_opt, x, V, L_opt, k_opt, z_opt);  % Applicazione vettoriale della funzione f

% Visualizzare il confronto tra i dati reali e quelli predetti
for i = 1:length(V)
    figure;  % Crea 5 grafici uno sopra l'altro
    % Interpolazione dei dati originali per il voltaggio V(i)
    F_interpolated = interp1(x, F_filtered(:,i), x, 'spline');
    plot(x, F_interpolated, 'o-', 'DisplayName', ['Dati reali (V = ' num2str(V(i)) ')']);
    hold on;
    plot(x, F_pred(:,i), 'x-', 'DisplayName', ['Forza predetta (V = ' num2str(V(i)) ')']);
    xlabel('Posizione x');
    ylabel('Forza F');
    legend;
    title(['Confronto forze per V = ', num2str(V(i))]);
end
