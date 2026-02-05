% PRINT DELLE BOBINE ATTIVE OGNI MM

lenght_mover = 150;
pitch_coils = 1;
dim_magnet = [10 15];
diam_coil = 25;
num_coils = 7;
lenght_tot = num_coils*diam_coil + (num_coils-1)*pitch_coils;
pitch = diam_coil + pitch_coils;
dx = 0.5;
dx1 = 0.1;
x_values = 0 : dx : 26*7;
test_values = 0 : dx1 : 23;

% part that interpolates the values and find an approximation of the force
% import variables d, attractionfx, repulsionfx, offfx

F_att = fit(d, F24VATT, 'cubicspline');
F_rep = fit(d, F24VREP, 'cubicspline');


Values_attraction = F_att(test_values);
Values_repultion = F_rep(test_values);


F_x = zeros(1, length(x_values)); 
num_magn = 4;

for k = 33 : 1 : 3  % Minimum pitch equal to coil pitch, max double
    for p = 72 : 1 :  72
        for f =  114 : 1 : 114
        
        % Apri il file per la scrittura
        fid = fopen('output.txt', 'w');

        % Scrivi l'intestazione
        fprintf(fid, 'X\tCOIL_ATTRACTION\tCOIL_REPULSION\n');
            
            
        mover = [0, k, p, f];
        coils = pitch * (0 : 1 : num_coils);
        o = 1;
        flag = 0;
        a = 0;
        scatti = zeros(5, length(x_values));
        total_forces = 0; % Initialize total force sum for averaging
            %simulating the interaction between the magnet and the coils and
            %store the force value
    
            %determine the distance of the magn -1 for repulsion and
            %+1 for attraction for only the magnets which are interacting
    
            %update the magn vector
        for x = x_values
            a = a + 1;
            distance = zeros(num_magn, num_coils); %matrix storing the distance of each magnet to each coil
            mode = zeros(1, num_coils); %array containig the mode of the coil
            for s = 1 : num_coils  
                for j = 1 : num_magn
                     distance(j, s) = coils(s) - mover(j);
                    
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
    
         coil_attraction = find(mode == 1);
         coil_repulsion = find(mode == -1);

    % Converti gli array in stringhe per la scrittura nel file
    coil_attr_str = strjoin(string(coil_attraction), ',');
    coil_rep_str = strjoin(string(coil_repulsion), ',');

    % Scrivi i dati nel file
    fprintf(fid, '%f\t%s\t%s\n', x, coil_attr_str, coil_rep_str);




for s = 1:num_coils
                F_total_rep = 0;
                F_total_att = 0;
                
                for j = 1:num_magn
                    if abs(distance(j, s)) < 23  % Considero solo le distanze note
                        if distance(j, s) < 0
                            F_total_rep = F_total_rep - F_rep(abs(distance(j, s))); % Repulsione
                            F_total_att = F_total_att - F_att(abs(distance(j, s))); % Repulsione
                        else
                            F_total_rep = F_total_rep + F_rep(abs(distance(j, s))); % Attrazione
                            F_total_att = F_total_att + F_att(abs(distance(j, s))); % Attrazione
                        end
                    end
                end
                
                % Seleziono la modalità che genera la forza maggiore in valore assoluto
                if F_total_att > F_total_rep
                    mode(s) = 1;  % Attrazione
                else
                    mode(s) = -1; % Repulsione
                end
            end



        
           
            % pause(0.001);
        
            F_x(o) = total_force;
            o = o + 1;
            mover = mover + dx;
        end
        
        % Chiudi il file
        fclose(fid);

        disp('File output.txt generato con successo!');

            
    end
    end
end