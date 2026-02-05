%programma che mostra l'andamento di accelerazone, decelerazione, forza di
%stop e forza non controllabile(dovuta a interazione tra magnete e ferro)
%in funzione di x, anche in questo caso importare i dati del file cubicPm
%per far funzionare il tutto

lenght_mover = 150;
pitch_coils = 1;
dim_magnet = [10 15];
diam_coil = 25;
num_coils = 20;
lenght_tot = num_coils*diam_coil + (num_coils-1)*pitch_coils;
pitch = diam_coil + pitch_coils;
dx = 0.1;
x_values = 0 : dx : 26;
m = 0.7; %mass of the cart


% part that interpolates the values and find an approximation of the force
% import variables d, attractionfx, repulsionfx, offfx

F_att = fit(d, F24VATT, 'cubicspline');
F_rep = fit(d, F24VREP, 'cubicspline');
F_off = fit(d, FxOff, 'cubicspline');
% Creazione del vettore di punti per il plot
d_plot = linspace(min(d), max(d), 1000);
 
F_x_dec = zeros(2, length(x_values)); 
F_x_stop = zeros(2, length(x_values)); 
num_magn = 4;
i = 1;
%calcolo della forze 
for k = 46 : 1 : 46  
    for p = 92 : 1 : 92
        for b = 138 : 1 : 138
            for v = 23 : 1 : 23
            [flag, F_tot_acc] = compute_force_acc(k, p, b, v, F_att, F_rep);
            if flag == 1
                continue;
            end
            [flag, F_tot_dec] = compute_force_dec(k, p, b, v, F_att, F_rep);
            if flag == 1
                continue;
            end
            [flag, F_tot_stop] = compute_force_stop(k, p, b, v, F_att, F_rep);
            if flag == 1
                continue;
            end

            [flag, F_tot_nc] = compute_force_non_controllable(k, p, b, v, F_att);
            

            figure;
            plot(x_values, F_tot_nc, 'k', 'LineWidth', 2);
            grid on;
            box on;
            xlabel('X-position');
            ylabel('Force applied to the system');
            title(['Cogging Force with pitch ', num2str(k), ' and ', num2str(p), ' and ', num2str(b), ' and offset ', num2str(v)]);
            ylim([-15, 15]);  % Fissa l'intervallo dell'asse Y
            hold on;



            figure;
            plot(x_values, F_tot_dec, 'k', 'LineWidth', 2);
            grid on;
            box on;
            xlabel('X-position');
            ylabel('Force applied to the system');
            title(['Force DECELERATION with pitch ', num2str(k), ' and ', num2str(p), ' and ', num2str(b), ' and offset ', num2str(v)]);
            hold on;

            figure;
            plot(x_values, F_tot_stop, 'k', 'LineWidth', 2);
            grid on;
            box on;
            xlabel('X-position');
            ylabel('Force applied to the system');
            title(['Force STOP with pitch ', num2str(k), ' and ', num2str(p), ' and ', num2str(b), ' and offset ', num2str(v)]);


            figure;
            plot(x_values, F_tot_acc, 'k', 'LineWidth', 2);
            grid on;
            box on;
            xlabel('X-position');
            ylabel('Force applied to the system');
            title(['Force ACCELERATION with pitch ', num2str(k), ' and ', num2str(p), ' and ', num2str(b), ' and offset ', num2str(v)]);
            hold on;
            fprintf('configurazione da salvare: magnete 2: %d // magnete 3: %d // magnete 4: %d // offset: %d \n', k, p, b, v);

            end
        end
    end
end
% %interpolo risultante
    % F_interp_stop = fit(x_values(:), F_tot_stop(:), 'smoothingspline', 'SmoothingParam', 0.97);
    % Values_interpolation3 = F_interp_stop(x_values);
    % 
    % 
    % %Plot approximation
    % figure;
    % plot(x_values, Values_interpolation3);
    % grid on;
    % box on;
    % xlabel('X-position');
    % ylabel('Force');
    % title('Force on the mover STOP approximation');
    % hold on;
    % plot(x_values, F_tot_stop);

% %interpolo risultante
    % F_interp_acc = fit(x_values(:), F_tot_acc(:), 'smoothingspline', 'SmoothingParam', 0.97);
    % Values_interpolation = F_interp_acc(x_values);
    % 
    % 
    % %Plot approximation
    % figure;
    % plot(x_values, Values_interpolation);
    % grid on;
    % box on;
    % xlabel('X-position');
    % ylabel('Force');
    % title('Force ACCELERATION approximation');
    % hold on;
    % plot(x_values, F_tot_acc);


% %interpolo risultante
    % F_interp_dec = fit(x_values(:), F_tot_dec(:), 'smoothingspline', 'SmoothingParam', 0.97);
    % Values_interpolation2 = F_interp_dec(x_values);
    % 
    % 
    % %Plot approximation
    % figure;
    % plot(x_values, Values_interpolation2);
    % grid on;
    % box on;
    % xlabel('X-position');
    % ylabel('Force');
    % title('Force DECELERATION approximation');
    % hold on;
    % plot(x_values, F_tot_dec);


% %DINAMICA DEL SISTEMA ATTRAVERSO CALCOLO NUMERICO
    % %STO CERCANDO V(t) NOTO CHE a = a(x(t)) ATTRAVERSO LA 2LDN
    % %E CHE x(t) = 0.5a(x(t))t^2 + v(t)*t  
    % %voglio che si fermi in un punto
    % 
    % 
    % F_1 = @(x) feval(F_interp_acc, x);
    % F_2 = @(x) feval(F_interp_dec, x);
    % T = 26;
    % F_acc = @(x) F_1(mod(x, T)); %forza periodica, di periodo spaziale 26mm
    % F_dec = @(x) F_2(mod(x, T));
    % 
    % tspan = linspace(0, 0.2, 1000); 
    % 
    % 
    % xdd = @(x) F_acc(x) / m;
    % dynamics = @(t, y) [y(2); xdd(y(1)*1e3)];
    % x_0 = 0;
    % v_0 = 0;
    % y0 = [x_0, v_0];
    % 
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



function [flag, F_tot_acc] = compute_force_acc(k, p, b, v, F_att, F_rep)
    pitch_coils = 1;
    diam_coil = 25;
    num_coils = 20;
    pitch = diam_coil + pitch_coils;
    dx = 0.1;
    x_values = 0 : dx : 26;
    F_x_acc = zeros(2, length(x_values));
    F_tot_acc = zeros(1, length(x_values));

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
                            else
                                Fx = F_rep(abs(distance(j, s)));
                            end
                        else
                            if distance(j, s) < 0
                                Fx = -F_att(abs(distance(j, s)));
                            else
                                Fx = F_att(abs(distance(j, s)));
                            end
                        end
                        total_force = total_force + Fx;
                    end
                end
            end          
            F_x_acc(i, o) = total_force;
        end
        temp = F_x_acc(1, o) + F_x_acc(2, o);
        % if temp < 5
        %     flag = 1;
        %     return;
        % end
        o = o + 1;
        mover = mover + dx;
    end
    F_tot_acc = F_x_acc(1, :) + F_x_acc(2, :);
end



function [flag, F_tot_dec] = compute_force_dec(k, p, b, v, F_att, F_rep)
    pitch_coils = 1;
    diam_coil = 25;
    num_coils = 20;
    pitch = diam_coil + pitch_coils;
    dx = 0.1;
    x_values = 0:dx:26;
    F_x_acc = zeros(2, length(x_values));
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
                            else
                                Fx = F_rep(abs(distance(j, s)));
                            end
                        else
                            if distance(j, s) < 0
                                Fx = -F_att(abs(distance(j, s)));
                            else
                                Fx = F_att(abs(distance(j, s)));
                            end
                        end
                        total_force = total_force + Fx;
                    end
                end
            end 
            F_x_acc(i, o) = total_force;
        end
        temp = F_x_acc(1, o) + F_x_acc(2, o);
        % if temp > 0
        % flag = 1;
        % return;
        % end
        o = o + 1;
        mover = mover + dx;
    end
    F_tot_dec = F_x_acc(1, :) + F_x_acc(2, :);
end




function [flag, F_tot_stop] = compute_force_stop(k, p, b, v, F_att, F_rep)
    pitch_coils = 1;
    diam_coil = 25;
    num_coils = 20;
    pitch = diam_coil + pitch_coils;
    dx = 0.1;
    x_values = 0:dx:26;
    F_x_acc = zeros(2, length(x_values));
    F_tot_stop = zeros(1, length(x_values));

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
                            else
                                Fx = F_rep(abs(distance(j, s)));
                            end
                        else
                            if distance(j, s) < 0
                                Fx = -F_att(abs(distance(j, s)));
                            else
                                Fx = F_att(abs(distance(j, s)));
                            end
                        end
                        total_force = total_force + Fx;
                    end
                end
            end
            F_x_acc(i, o) = total_force;
         end
        temp = abs(F_x_acc(1, o) + F_x_acc(2, o));
        % if temp > 4.5
        %     flag = 1;
        %     return;
        % end
        o = o + 1;
        mover = mover + dx;
    end
    F_tot_stop = F_x_acc(1, :) + F_x_acc(2, :);
end



function [flag, F_tot_nc] = compute_force_non_controllable(k, p, b, v, F_off)
    pitch_coils = 1;
    diam_coil = 25;
    num_coils = 20;
    pitch = diam_coil + pitch_coils;
    dx = 0.1;
    x_values = 0 : dx : 26;
    F_x_acc = zeros(2, length(x_values));
    F_tot_nc = zeros(1, length(x_values));

    mover = [0, k, p, b; v, k + v, p + v, b + v];    
    coils = pitch * (0:num_coils);
    o = 1;
    flag = 1;
    
    for x = x_values
        for i = 1:2
            total_force = 0;
            for j = 1:4
                for s = 1:num_coils
                    distance = coils(s) - mover(i, j);
                    Fx = 0; % Inizializza Fx per evitare problemi
                    if abs(distance) < 23
                       if distance < 0
                             Fx = -F_off(abs(distance));
                       else
                             Fx = F_off(abs(distance));
                       end
                    end
                    total_force = total_force + Fx;
                end
            end          
            F_x_acc(i, o) = total_force;
        end
        o = o + 1;
        mover = mover + dx;
    end
    F_tot_nc = F_x_acc(1, :) + F_x_acc(2, :);
end
