%programma che calcola per quanto tempo devono stare accesi/repellere i
%primi 7 solenoidi, usato per preparare arduino, 

lenght_mover = 150;
pitch_coils = 1;
dim_magnet = [10 15];
diam_coil = 25;
num_coils = 20;
lenght_tot = num_coils*diam_coil + (num_coils-1)*pitch_coils;
pitch = diam_coil + pitch_coils;
dx = 0.1;
x_values = 0 : dx : 26;
test_values = 0 : dx : 20;
m = 0.778; %mass of the cart


% part that interpolates the values and find an approximation of the force
% import variables d, attractionfx, repulsionfx, offfx

F_att = fit(d, FxAttractioncut, 'cubicspline');
F_rep = fit(d, FxRepultioncut, 'cubicspline');


Values_attraction = F_att(test_values);
Values_repultion = F_rep(test_values);



F_x = zeros(2, length(x_values)); 
num_magn = 4;

%calcolo della forza sul mover

for k = 33 : 1 : 33  
    for p = 72 : 1 : 72
        for b = 117 : 1 : 117
        for v = 3 : 1 : 3
        mover = zeros(2, 4);  
        mover(1, :) = [0, k , p, b];
        mover(2, :) = [v, k + v, p + v, b + v];
        coils = pitch * (0 : 1 : num_coils);
        fav = 0;
        flag = 0;
        for i = 1 : 1 : 2  
        o = 1;
        flag = 0;
        for x = x_values
            distance = zeros(num_magn, num_coils); %matrix storing the distance of each magnet to each coil
            mode = zeros(1, num_coils); %array containig the mode of the coil
            for s = 1 : num_coils  
                for j = 1 : num_magn
                    distance(j, s) = coils(s) - mover(i, j);
                 end
            end
            %the matrix contains the distance and its sign, if dist <0
            %repulsion, if dist >0 attraction
            % if abs (distance) > 20 no influence on the system
                
            % now I want to evaluate the force between magnets and coils I have to look inside the matrix
            % and understand the forces on a single magnet, then sum the forces of all the magnets
            sum = 0;
            total_force = 0;
            
                       %define the mode of the coils for each iteration
            
            for s = 1:num_coils
                mode_coil = 0;
                for j = 1:num_magn
                    if abs(distance(j, s)) < 20
                    mode_coil = sign(distance(j, s)); 
                    if mode_coil > 0
                        if j < num_magn && abs(distance(j + 1, s)) < 20
                        % Select the value with the smaller magnitude, keeping its sign
                        mode_coil = sign(distance(j + (abs(distance(j + 1, s)) < abs(distance(j, s))), s));
                        break;
                        end
                    end
                    if mode_coil < 0
                        if abs(distance(j, s)) < 15
                        if j < num_magn && abs(distance(j + 1, s)) < 20
                        % Select the value with the smaller magnitude, keeping its sign
                        mode_coil = sign(distance(j + (abs(distance(j + 1, s)) < abs(distance(j, s))), s));
                        break;
                        end
                        else
                            mode_coil = 0;
                        end
                    end
                    end
                end
                
                mode(1, s) = mode_coil;
            end

            for j = 1 : num_magn
                for s = 1 : num_coils
                    if abs(distance(j, s)) < 20     %no information if bigger than 20
                        if mode(s) < 0 
                            Fx = - F_rep(abs(distance(j, s))); %repultion negative
                        else 
                            Fx = F_att(abs(distance(j, s))); %attraction positive
                        end
                        total_force = total_force + Fx;
                    end
                end
            end
            if total_force < 0
               flag = 1;
               break;
            end

        F_x(i, o) = total_force;
        o = o + 1;
        mover(i, :) = mover(i, :) + dx;
        
        end
        end

        F_tot = F_x(1, :) + F_x(2, :);
        minValue = min(F_tot); 
        
        if minValue < 8.5
            flag = 1;
        end


        end
        end
    end
end


%interpolo risultante
F_interp = fit(x_values(:), F_tot(:), 'smoothingspline', 'SmoothingParam', 0.97);
Values_interpolation = F_interp(x_values);

filename = 'output.txt'; % Assicurati che il file sia nella cartella di lavoro
[solenoid_matrix, X_values] = process_solenoids(filename);


%DINAMICA DEL SISTEMA ATTRAVERSO CALCOLO NUMERICO
%STO CERCANDO V(t) NOTO CHE a = a(x(t)) ATTRAVERSO LA 2LDN
%E CHE x(t) = 0.5a(x(t))t^2 + v(t)*t  

F_1 = @(x) feval(F_interp, x);
T = 26;
F = @(x) F_1(mod(x, T)); %forza periodica, di periodo spaziale 26mm
xdd = @(x) F(x) / m;
dynamics = @(t, y) [y(2); xdd(y(1)*1e3)];
x_0 = X_values(1)*1e-3;
v_0 = 0;
y0 = [x_0, v_0];
tspan = linspace(0, 0.1, 10000); 
% Risoluzione numerica con ode45
[t, Y] = ode45(dynamics, tspan, y0);

% Estrazione delle soluzioni
x_t = Y(:, 1)*1e3;
v_t = Y(:, 2);
a_t = arrayfun(@(x) xdd(x), x_t*1e3); 

fid = fopen('tempi.txt', 'w'); % Apri il file per scrivere

[num_rows, num_coils] = size(solenoid_matrix); % Ottieni la dimensione della matrice

for j = 1:num_coils
    q = 1; % Inizializza l'indice di riga
    
    while q <= num_rows
        if solenoid_matrix(q, j) == 1
            start = X_values(q);
            while q < num_rows && solenoid_matrix(q, j) == 1
                q = q + 1;
            end
            finish = X_values(q);
           
            % Trova le posizioni di start e finish in x_t con il valore più vicino
            [~, start_idx] = min(abs(x_t - start));
            [~, finish_idx] = min(abs(x_t - finish));

            
            if ~isempty(start_idx) && ~isempty(finish_idx) % Se troviamo entrambi i valori
                start_time = t(start_idx);
                finish_time = t(finish_idx);
                
                % Scrivi sul file
                fprintf(fid, 'Solenoide %d in modalità ON da t = %.6f a t = %.6f\n', j, start_time, finish_time);
                
            end
        end

        if solenoid_matrix(q, j) == -1
            start = X_values(q);
            while q < num_rows && solenoid_matrix(q + 1, j) == -1
                q = q + 1;
            end
            finish = X_values(q);
            
            % Trova le posizioni di start e finish in x_t con tolleranza 1e-3
            % Trova le posizioni di start e finish in x_t con il valore più vicino
            [~, start_idx] = min(abs(x_t - start));
            [~, finish_idx] = min(abs(x_t - finish));

            
            if ~isempty(start_idx) && ~isempty(finish_idx) % Se troviamo entrambi i valori
                start_time = t(start_idx);
                finish_time = t(finish_idx);
                
                % Scrivi sul file
                fprintf(fid, 'Solenoide %d in modalità REP da t = %.6f a t = %.6f\n', j, start_time, finish_time);
                
            end
        end
        q = q + 1; % Passa alla prossima riga
    end
end

fclose(fid); % Chiudi il file

disp('File tempi.txt generato con successo!');



function [solenoid_matrix, X_values] = process_solenoids(filename)
    % Legge il file di testo
    data = readtable(filename, 'Delimiter', '\t', 'ReadVariableNames', true);
    
    % Estrae i valori di X in un vettore colonna
    X_values = data.X;
    
    % Trova tutti i solenoidi unici presenti
    attraction_cells = cellfun(@(x) str2num(x), data.COIL_ATTRACTION, 'UniformOutput', false);
    repulsion_cells = cellfun(@(x) str2num(x), data.COIL_REPULSION, 'UniformOutput', false);
    all_solenoids = unique([horzcat(attraction_cells{:}), horzcat(repulsion_cells{:})]);
    
    % Numero di righe del file
    num_rows = height(data);
    num_solenoids = max(all_solenoids); % Numero massimo di solenoidi identificati
    
    % Inizializza la matrice output
    solenoid_matrix = zeros(num_rows, num_solenoids);
    
    % Popola la matrice con 1 per attraction e -1 per repulsion
    for i = 1:num_rows
        if ~isempty(attraction_cells{i})
            solenoid_matrix(i, attraction_cells{i}) = 1;
        end
        if ~isempty(repulsion_cells{i})
            solenoid_matrix(i, repulsion_cells{i}) = -1;
        end
    end
end