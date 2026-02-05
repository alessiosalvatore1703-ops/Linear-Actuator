% Caricamento dati da Excel
filename = 'cubic_pm.xlsx';
data = readmatrix(filename);
x = data(2:end, 1);  % Posizioni x (escludendo l'intestazione)
V = [0, -12, 12, -24, 24];  % Valori di V (esclusa l'intestazione)
%F = data(2:end, 2:end);  % Forze (matrice con le forze per i vari V)

% Grado del polinomio di approssimazione (puoi modificarlo)
degree = 5;

% Creazione di 5 figure, una per ogni valore di V
for i = 1:length(V)
    % Estrai i dati di forza per il corrente V
    y = F(:, i);
    
    % Calcola i coefficienti del polinomio
    p = polyfit(x, y, degree);
    
    % Crea punti più fitti per una curva più liscia
    x_fine = linspace(min(x), max(x), 100);
    y_fine = polyval(p, x_fine);
    
    % Crea una nuova figura
    figure;
    
    % Plot dei dati originali
    plot(x, y, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on;
    
    % Plot dell'approssimazione polinomiale
    plot(x_fine, y_fine, 'r-', 'LineWidth', 2);
    
    % Aggiunta di titolo e legenda
    title(sprintf('Approssimazione polinomiale per V = %d', V(i)));
    xlabel('Posizione x');
    ylabel('Forza F');
    legend('Dati reali', 'Approssimazione polinomiale', 'Location', 'best');
    
    % Griglia per migliore leggibilità
    grid on;
    
    % Personalizzazione dell'aspetto
    set(gca, 'FontSize', 12);
end