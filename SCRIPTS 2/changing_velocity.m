%% Controllo del mover

lenght_mover = 150;
pitch_coils = 1;
dim_magnet = [10 15];
diam_coil = 25;
num_coils = 20;
lenght_tot = num_coils * diam_coil + (num_coils - 1) * pitch_coils;
pitch = diam_coil + pitch_coils;
dx = 0.1;
x_values = 0 : dx : 26;
test_values = 0 : dx : 20;
m = 1; % Massa del carrello

% Interpolazione dei dati
F_att = fit(d, FxAttractioncut, 'cubicspline');
F_rep = fit(d, FxRepultioncut, 'cubicspline');
Values_attraction = feval(F_att, test_values);
Values_repulsion = feval(F_rep, test_values);
F_x = zeros(2, length(x_values)); 
num_magn = 4;

% Determinazione del coefficiente k tramite regressione lineare
p = polyfit(V, FzRep, 1); % Fit lineare (y = p(1)*x + p(2))
k = p(1);
V_values = 0 : dx : 24;
p_values = polyval(p, V_values);




% Calcolo della forza media F_off
F_off = zeros(size(test_values));
for i = 1:length(test_values)
    F_off(i) = (Values_attraction(i) + Values_repulsion(i)) / 2;
end

% Definizione della funzione F_V
F_V = @(V) F_off + k * V;

% Plot della forza per vari valori di V
figure;
grid on;
colormap = jet(length(-24:4:24)); % Creazione di una mappa di colori
index = 1;
for j = -24 : 4 : 24
    hold on;
    plot(test_values, F_V(j), 'Color', colormap(index, :));
    index = index + 1;
end
plot(test_values, Values_attraction, 'g');
plot(test_values, Values_repulsion, 'b');
hold off;

figure;
grid on;
plot(V, FzRep);
hold on;
plot(V_values, p_values);

% %%
% %%calcolo della forza sul mover
% 
% for k = 33 : 1 : 33  
%     for p = 72 : 1 : 72
%         for b = 117 : 1 : 117
%         for v = 3 : 1 : 3
%         mover = zeros(2, 4);  
%         mover(1, :) = [0, k , p, b];
%         mover(2, :) = [v, k + v, p + v, b + v];
%         coils = pitch * (0 : 1 : num_coils);
%         fav = 0;
%         flag = 0;
%         %a = 0;
%         %scatti = zeros(5, length(x_values));
%             %simulating the interaction between the magnet and the coils and
%             %store the force value
% 
%             %determine the distance of the magn -1 for repulsion and
%             %+1 for attraction for only the magnets which are interacting
%             %update the magn vector
%         for i = 1 : 1 : 2  
%         o = 1;
%         flag = 0;
%         for x = x_values
%             %a = a + 1;
%             distance = zeros(num_magn, num_coils); %matrix storing the distance of each magnet to each coil
%             mode = zeros(1, num_coils); %array containig the mode of the coil
%             for s = 1 : num_coils  
%                 for j = 1 : num_magn
%                     distance(j, s) = coils(s) - mover(i, j);
%                  end
%             end
%             %the matrix contains the distance and its sign, if dist <0
%             %repulsion, if dist >0 attraction
%             % if abs (distance) > 20 no influence on the system
% 
%             % now I want to evaluate the force between magnets and coils I have to look inside the matrix
%             % and understand the forces on a single magnet, then sum the forces of all the magnets
%             sum = 0;
%             total_force = 0;
% 
%             %define the mode of the coils for each iteration
% 
%             for s = 1:num_coils
%                 mode_coil = 0;
%                 for j = 1:num_magn
%                     if abs(distance(j, s)) < 20
%                     mode_coil = sign(distance(j, s)); 
%                         if j < num_magn && abs(distance(j + 1, s)) < 20
%                         % Select the value with the smaller magnitude, keeping its sign
%                         mode_coil = sign(distance(j + (abs(distance(j + 1, s)) < abs(distance(j, s))), s));
%                         break;
%                         end
%                     end
%                 end
%                 mode(1, s) = mode_coil;
%             end
% 
% 
% 
%             for j = 1 : num_magn
%                 for s = 1 : num_coils
%                     if abs(distance(j, s)) < 20     %no information if bigger than 20
%                         if mode(s) < 0 
%                             Fx = - F_rep(abs(distance(j, s))); %repultion negative
%                         else 
%                             Fx = F_att(abs(distance(j, s))); %attraction positive
%                         end
%                         total_force = total_force + Fx;
%                     end
%                 end
%             end
%             if total_force < 0
%                flag = 1;
%                break;
%             end
% 
%        %plot magneti 
%         %     clf;
%         %     hold on;
%         %     axis([0, lenght_tot, -10, 10]);
%         %     title(['Magnets with pitch ', num2str(k)])
%         % 
%         %     display_radius = 25; % Raggio dei punti
%         % 
%         % % Disegna solenoidi
%         % for s = 1:num_coils
%         %     color = [1 1 1]; % Bianco di default
%         %     if mode(1, s) == -1
%         %         color = [0 1 0]; % Blu per repulsione
%         %     end
%         %     if mode(1, s) == 1
%         %         color = [1 0 0]; % Rosso per attrazione
%         %     end
%         %     scatter(coils(s), 0, 100, color, 'filled'); % Solenoidi come punti
%         % end
%         % 
%         % % Disegna magneti
%         % for j = 1:num_magn
%         %     scatter(mover(j), 0, 100, [0 0 0], 'filled'); % Magneti come punti neri
%         % end
%         % 
%         % pause(0.001);
% 
%         F_x(i, o) = total_force;
%         o = o + 1;
%         mover(i, :) = mover(i, :) + dx;
% 
%         end
%         end
% 
%         F_tot = F_x(1, :) + F_x(2, :);
%         minValue = min(F_tot); 
% 
%         if minValue < 8.5
%             flag = 1;
%         end
% 
%         if flag == 0
%         % Plot
%         figure;
%         plot(x_values, F_tot, 'k', 'LineWidth', 2); % Forza in nero più visibile
%         grid on;
%         box on;
%         xlabel('X-position');
%         ylabel('Force applied to the system');
%         title(['Force vs displacement with pitch ', num2str(k), ' and ', num2str(p),' and ', num2str(b) ' e disassamento', num2str(v), 'valore minimo' num2str(minValue)]);
%         hold on;
% 
%         plot(x_values, F_x(1, :), "r");
%         plot(x_values, F_x(2, :), "g");
% 
%         end
% 
% 
% 
%         end
%             end
%         end
%     end
% 
% 
% %interpolo risultante
% F_interp = fit(x_values(:), F_tot(:), 'smoothingspline', 'SmoothingParam', 0.97);
% Values_interpolation = F_interp(x_values);
% 
% 
% %Plot approximation
% figure;
% plot(x_values, Values_interpolation);
% grid on;
% box on;
% xlabel('X-position');
% ylabel('Force');
% title('Force on the mover approximation');
% hold on;
% plot(x_values, F_tot);
% 
% %DINAMICA DEL SISTEMA ATTRAVERSO CALCOLO NUMERICO
% %STO CERCANDO V(t) NOTO CHE a = a(x(t)) ATTRAVERSO LA 2LDN
% %E CHE x(t) = 0.5a(x(t))t^2 + v(t)*t  
% 
% F_1 = @(x) feval(F_interp, x);
% T = 26;
% F = @(x) F_1(mod(x, T)); %forza periodica, di periodo spaziale 26mm
% xdd = @(x) F(x) / m;
% dynamics = @(t, y) [y(2); xdd(y(1)*1e3)];
% x_0 = 3;
% v_0 = 0;
% y0 = [x_0, v_0];
% tspan = linspace(0, 1, 1000); 
% % Risoluzione numerica con ode45
% [t, Y] = ode45(dynamics, tspan, y0);
% 
% % Estrazione delle soluzioni
% x_t = Y(:, 1);
% v_t = Y(:, 2);
% a_t = arrayfun(@(x) xdd(x), x_t*1e3); 
% 
% % Plot dei risultati
% figure;
% 
% % Subplot della posizione
% subplot(3,1,1);
% plot(t, x_t, 'b', 'LineWidth', 1.5);
% xlabel('Tempo [s]');
% ylabel('Posizione x(t)');
% title('Evoluzione della posizione');
% grid on;
% 
% % Subplot della velocità
% subplot(3,1,2);
% plot(t, v_t, 'r', 'LineWidth', 1.5);
% xlabel('Tempo [s]');
% ylabel('Velocità v(t) ');
% title('Evoluzione della velocità');
% grid on;
% 
% % Subplot dell'accelerazione
% subplot(3,1,3);
% plot(t, a_t, 'g', 'LineWidth', 1.5);
% xlabel('Tempo [s]');
% ylabel('Accelerazione a(t)');
% title('Evoluzione accelerazione');
% grid on;