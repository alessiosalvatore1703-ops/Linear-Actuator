%programma che calcola passo ottimale 3 magneti, passi non regolari
%considero mover largo 100mm e lungo 150mm
%passo solenoidi 1mm

lenght_mover = 150;
pitch_coils = 1;
dim_magnet = [10 15];
diam_coil = 25;
num_coils = 10;
lenght_tot = num_coils*diam_coil + (num_coils-1)*pitch_coils;
pitch = diam_coil + pitch_coils;
dx = 0.01;
x_values = 0 : dx : 26;
test_values = 0 : dx : 23;

% part that interpolates the values and find an approximation of the force
% import variables d, attractionfx, repulsionfx, offfx

F_att = fit(d, F24VATT, 'cubicspline');
F_rep = fit(d, F24VREP, 'cubicspline');


Values_attraction = F_att(test_values);
Values_repultion = F_rep(test_values);


%Plot approximations
% figure;
% plot(test_values, Values_attraction);
% grid on;
% box on;
% xlabel('X-position');
% ylabel('Force');
% title('Force vs displacement polynomial approximation Attraction');
% hold on;
% plot(d, F_att1);
% 
% figure;
% plot(test_values, Values_repultion);
% grid on;
% box on;
% xlabel('X-position');
% ylabel('Force');
% title('Force vs displacement polynomial approximation Repulsion');
% hold on;
% plot(d, F_rep1);


F_x = zeros(1, length(x_values)); 
num_magn = 4;

for k = 39 : 1 : 39  % Minimum pitch equal to coil pitch, max double
    for p = 72 : 1 : 72 
        for f =  110 : 1 : 110
            
        mover = [0, k , p, f];
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
                    
                    if (abs(distance(j, s) - 5) < 1e-2 && distance(j, s) > 0) % massimo 1
                        scatti (1, a) = a;
                    end
                    if (abs(distance(j, s) - 8) < 1e-2 && distance(j, s) > 0) % minimo 1
                        scatti (2, a) = a;
                    end
                   
                    if (abs(abs(distance(j, s)) - 10) < 1e-2 && distance(j, s) < 0) % massimo 2
                        scatti (4, a) = a;
                    end
                    if (abs(abs(distance(j, s)) - 15) < 1e-2 && distance(j, s) < 0) % shell
                            scatti (5, a) = a;
                    end
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
 
            

            for j = 1 : num_magn
                for s = 1 : num_coils
                    if abs(distance(j, s)) < 23     %no information if bigger than 20
                        if mode(s) < 0
                            if distance(j, s) < 0
                            Fx = - F_rep(abs(distance(j, s))); % Repulsione
                            else
                            Fx = F_rep(abs(distance(j, s))); % Attrazione
                            end
                        end
                        if mode(s) > 0
                            if distance(j, s) < 0
                            Fx = - F_att(abs(distance(j, s))); % Repulsione
                            else
                            Fx = F_att(abs(distance(j, s))); % Attrazione
                            end
                        end

                        total_force = total_force + Fx;
                    end
                end
            end

            % if total_force < 2.5
            %     flag = 1;
            %     break;
            % end



        
            %     clf;
            %     hold on;
            %     axis([0, lenght_tot, -10, 10]);
            %     title(['Magnets with pitch ', num2str(k)])
            % 
            %     display_radius = 25; % Raggio dei punti
            % 
            % % Disegna solenoidi
            % for s = 1:num_coils
            %     color = [1 1 1]; % Bianco di default
            %     if mode(1, s) == -1
            %         color = [0 1 0]; % verde per repulsione
            %     end
            %     if mode(1, s) == 1
            %         color = [1 0 0]; % Rosso per attrazione
            %     end
            %     scatter(coils(s), 0, 100, color, 'filled'); % Solenoidi come punti
            % end
            % 
            % % Disegna magneti
            % for j = 1:num_magn
            %     scatter(mover(j), 0, 100, [0 0 0], 'filled'); % Magneti come punti neri
            % end
            % 
            % pause(0.001);
        
            F_x(o) = total_force;
            o = o + 1;
            mover = mover + dx;
        end
        

            if flag == 0
            % Plot
            figure;
            plot(x_values, F_x, 'k', 'LineWidth', 1.5); % Forza in nero più visibile
            grid on;
            box on;
            xlabel('X-position');
            ylabel('Force applied to the system');
            title(['Force vs displacement with pitch ', num2str(k), ' and ', num2str(p), ' and ', num2str(f)]);
            
            hold on;
            for w = 1:length(x_values)
                if scatti(1, w) ~= 0 
                    xline(scatti(1, w) * dx - dx, 'Color', [0 0 0.55], 'LineWidth', 3); %Guscio sx: Blu ([0 0 0.55])
                end
                if scatti(2, w) ~= 0 
                    xline(scatti(2, w) * dx - dx, 'Color', [0.5 0 0.5], 'LineWidth', 3); % Picco sx: Viola ([0.5 0 0.5])
                end 
                if scatti(3, w) ~= 0 
                    xline(scatti(3, w) * dx - dx, 'Color', [0.55 0 0], 'LineWidth', 3); % Nucleo centrale: Rosso ([0.55 0 0])
                end 
                if scatti(4, w) ~= 0 
                    xline(scatti(4, w) * dx - dx, 'Color', [0 1 0], 'LineWidth', 3); % Picco dx: Verde ([0 1 0])
                end 
                if scatti(5, w) ~= 0 
                    xline(scatti(5, w) * dx - dx, 'Color', [1 0.65 0], 'LineWidth', 3); %Guscio dx: Arancione ([1 0.65 0])
                end
            end 
        end
    end
    end
end