% PRINT DELLE BOBINE ATTIVE OGNI MM

lenght_mover = 150;
pitch_coils = 1;
dim_magnet = [10 15];
diam_coil = 25;
num_coils = 10;
lenght_tot = num_coils * diam_coil + (num_coils - 1) * pitch_coils;
pitch = diam_coil + pitch_coils;
dx = 0.5;
dx1 = 0.1;
x_values = 0 : dx : 38;
test_values = 0 : dx1 : 20;

% Import variables d, attractionfx, repulsionfx, offfx
F_att = fit(d, F_att1, 'cubicspline');
F_rep = fit(d, F_rep1, 'cubicspline');

Values_attraction = F_att(test_values);
Values_repulsion = F_rep(test_values);

F_x = zeros(1, length(x_values));
num_magn = 4;

% Variabili per il salvataggio dati
data.X = [];
data.COIL_ATTRACTION = {};
data.COIL_REPULSION = {};

for k = 30 : 1 : 30  % Minimum pitch equal to coil pitch, max double
    for p = 72 : 1 : 72
        for f = 114 : 1 : 114
            
            mover = [0, k, p, f];
            coils = pitch * (0 : 1 : num_coils);
            o = 1;
            a = 0;
            total_forces = 0; % Initialize total force sum for averaging
            
            for x = x_values
                a = a + 1;
                distance = zeros(num_magn, num_coils); % Matrice delle distanze magnete-bobina
                mode = zeros(1, num_coils); % Modalità della bobina
                
                for s = 1 : num_coils  
                    for j = 1 : num_magn
                        distance(j, s) = coils(s) - mover(j);
                    end
                end

                % Determina la modalità delle bobine (attrazione/repulsione)
                for s = 1:num_coils
                    mode_coil = 0;
                    for j = 1:num_magn
                        if abs(distance(j, s)) < 20
                            mode_coil = sign(distance(j, s));
                            if mode_coil > 0
                                if j < num_magn && abs(distance(j + 1, s)) < 20
                                    mode_coil = sign(distance(j + (abs(distance(j + 1, s)) < abs(distance(j, s))), s));
                                    break;
                                end
                            end
                            if mode_coil < 0
                                if abs(distance(j, s)) < 15
                                    if j < num_magn && abs(distance(j + 1, s)) < 20
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

                coil_attraction = find(mode == 1);
                coil_repulsion = find(mode == -1);

                % Salvataggio nei dati
                data.X = [data.X, x];
                data.COIL_ATTRACTION{end+1} = coil_attraction;
                data.COIL_REPULSION{end+1} = coil_repulsion;

                F_x(o) = total_forces;
                o = o + 1;
                mover = mover + dx;
            end
        end
    end
end

N = length(data.X);
num_coils = 10;  % Assicurati che questo valore sia coerente con il tuo modello

% Prealloca le matrici (dimensione: [tempo x bobina])
attraction_matrix = zeros(N, num_coils);
repulsion_matrix  = zeros(N, num_coils);

for idx = 1:N
    % Legge i vettori di bobine attive per ciascun tempo
    coilAttr = data.COIL_ATTRACTION{idx};
    coilRep  = data.COIL_REPULSION{idx};
    
    % Imposta a 1 la posizione corrispondente alla bobina se è attiva
    if ~isempty(coilAttr)
        attraction_matrix(idx, coilAttr) = 1;
    end
    if ~isempty(coilRep)
        repulsion_matrix(idx, coilRep) = 1;
    end
end

% CREAZIONE DEI TIMESERIES
ts_attraction = timeseries(attraction_matrix, data.X);
ts_repulsion  = timeseries(repulsion_matrix, data.X);

% Aggiorna la struttura dati: sostituisci le celle con le timeseries
data.COIL_ATTRACTION = ts_attraction;
data.COIL_REPULSION  = ts_repulsion;

% Salva il file .mat con i dati formattati per Simulink
save('output.mat', 'data');
disp('File output.mat generato con successo!');