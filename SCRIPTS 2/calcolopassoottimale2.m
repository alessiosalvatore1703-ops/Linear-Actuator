%programma che calcola passo ottimale
%considero mover largo 100mm e lungo 150mm
%passo solenoidi 1mm

lenght_mover = 150;
pitch_coils = 1;
dim_magnet = [10 15];
diam_coil = 25;
num_coils = 10;
lenght_tot = num_coils*diam_coil + (num_coils-1)*pitch_coils;
pitch = diam_coil + pitch_coils;
dx = 0.1;
x_values = 0 : dx : 26;
test_values = 0 : dx : 20;

% part that interpolates the values and find an approximation of the force
% import variables d, attractionfx, repulsionfx, offfx

F_att = fit(d, FxAttractioncut, 'cubicspline');
F_rep = fit(d, FxRepultioncut, 'cubicspline');


Values_attraction = F_att(test_values);
Values_repultion = F_rep(test_values);


% Plot approximations
% figure;
% plot(test_values, Values_attraction);
% grid on;
% box on;
% xlabel('X-position');
% ylabel('Force');
% title('Force vs displacement polynomial approximation Attraction');
% hold on;
% plot(d, FxAttractioncut);
% 
% figure;
% plot(test_values, Values_repultion);
% grid on;
% box on;
% xlabel('X-position');
% ylabel('Force');
% title('Force vs displacement polynomial approximation Repulsion');
% hold on;
% plot(d, FxRepultioncut);


F_x = zeros(1, length(x_values));  

for k = 36 : 1 : 45  % Minimum pitch equal to coil pitch, max double
    num_magnets = ((lenght_mover - dim_magnet(1, 2)) / k) + 1; % look for the maximum amount of magnets that fits
    num_magn = ceil(num_magnets) - 1; %round
    gap = lenght_mover - (num_magn - 1) * k - dim_magnet(1, 2); % check that this displacement fits on the mover

    if gap > 0
        mover = 0 : k : (num_magn - 1) * k;
        coils = pitch * (0 : 1 : num_coils);
        o = 1;
        flag = 0;
        a = 0;
        scatti = zeros(5, length(x_values));
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
                    if (abs(distance(j, s) - 11.8) < 1e-2 && distance(j, s) > 0) % external shell left
                        scatti (1, a) = a;
                    end
                    if (abs(distance(j, s) - 2.5) < 1e-2 && distance(j, s) > 0) % peak attraction
                        scatti (2, a) = a;
                    end
                    if (abs(distance(j, s)) < 1e-2)  %nucleo centrale
                            scatti (3, a) = a; 
                    end
                    if (abs(abs(distance(j, s)) - 1.9) < 1e-2 && distance(j, s) < 0) % peak repulsion
                        scatti (4, a) = a;
                    end
                    if (abs(abs(distance(j, s)) - 11.4) < 1e-2 && distance(j, s) < 0) % external shell right
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
                mode_coil = 0;
                for j = 1:num_magn
                    if abs(distance(j, s)) < 20
                    mode_coil = sign(distance(j, s)); 
                        if j < num_magn && abs(distance(j + 1, s)) < 20
                        % Select the value with the smaller magnitude, keeping its sign
                        mode_coil = sign(distance(j + (abs(distance(j + 1, s)) < abs(distance(j, s))), s));
                        break;
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
        %         color = [0 1 0]; % Blu per repulsione
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
        title(['Force vs displacement with pitch ', num2str(k), ' and with ', num2str(num_magn), ' magnets']);
        
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
