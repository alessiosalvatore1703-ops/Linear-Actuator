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
test_values = 0 : dx : 23;
% part that interpolates the values and find an approximation of the force
% import variables d, attractionfx, repulsionfx, offfx

%Imposta un livello di smoothing più alto per avere curve più lisce
smoothing_param = 0.85; % Regola questo valore per più o meno smoothing

options = fitoptions('smoothingspline', 'SmoothingParam', smoothing_param);

% Applica la smoothing spline con il parametro scelto
F_att = fit(d, F_att1, 'cubicspline');
F_rep = fit(d, F_rep1, 'cubicspline');
F_off = fit(d, F_off1, 'cubicspline');

Values_attraction = F_att(test_values);
Values_repultion = F_rep(test_values);
Values_off = F_off(test_values);

% Plot approximations
figure;
plot(test_values, Values_attraction);
grid on;
box on;
xlabel('X-position');
ylabel('Force');
title('Force vs displacement polynomial approximation Attraction');
hold on;
plot(d, F_att1);

figure;
plot(test_values, Values_repultion);
grid on;
box on;
xlabel('X-position');
ylabel('Force');
title('Force vs displacement polynomial approximation Repulsion');
hold on;
plot(d, F_rep1);

figure;
plot(test_values, Values_off);
grid on;
box on;
xlabel('X-position');
ylabel('Force');
title('Force vs displacement polynomial approximation OFF');
hold on;
plot(d, F_off1);


F_x = zeros(1, length(x_values));  

for k = 10 : 1 : 52  % Minimum pitch equal to coil pitch, max double
    num_magnets = ((lenght_mover - dim_magnet(1, 2)) / k) + 1; % look for the maximum amount of magnets that fits
    num_magn = ceil(num_magnets) - 1; %round
    gap = lenght_mover - (num_magn - 1) * k - dim_magnet(1, 1); % check that this displacement fits on the mover

    if gap > 0
        mover = 0 : k : (num_magn - 1) * k;
        coils = pitch * (0 : 1 : num_coils);
        o = 1;
        flag = 0;
            %simulating the interaction between the magnet and the coils and
            %store the force value
    
            %determine the distance of the magn -1 for repulsion and
            %+1 for attraction for only the magnets which are interacting
    
            %update the magn vector
        for x = x_values
            distance = zeros(num_magn, num_coils); %matrix storing the distance of each magnet to each coil
            % Definizione delle interazioni tra solenoidi e magneti
            interaction = zeros(num_coils, num_magn);
            mode = zeros(1, num_coils); %array containig the mode of the coil
            for s = 1 : num_coils  
                for j = 1 : num_magn
                    distance(j, s) = coils(s) - mover(j);
                end
            end
            sum = 0;
            total_force = 0;
            %cerco config ideale
            F_total_rep = 0;
            F_total_att = 0;
            for s = 1 : num_coils
                for n = 1 : 2
                    Fb = 0;
                    Fc = 0;
                    for j = 1 : num_magn
                        if n == 1
                        if abs(distance(j, s)) < 23     %no information if bigger than 23
                            if distance(j, s) < 0
                                Fb = - F_rep(abs(distance(j, s))); %want to be attracted
                            end
                            if distance(j, s) > 0
                                Fb = F_rep(abs(distance(j, s))); %want to be rejected
                            end
                        F_total_rep = F_total_rep + Fb;
                        end
                        end
                        if n==2
                        if abs(distance(j, s)) < 23     %no information if bigger than 23
                            if distance(j, s) < 0
                                Fc = F_att(abs(distance(j, s))); %want to be attracted
                            end
                            if distance(j, s) > 0
                                Fc = - F_att(abs(distance(j, s))); %want to be rejected
                            end
                        F_total_att = F_total_att + Fc;
                        end                            
                        end
                    end
                end

                mode(s) = sign(F_total_att - F_total_rep); %
             end 
            



            for j = 1 : num_magn
                for s = 1 : num_coils
                    if abs(distance(j, s)) < 23     %no information if bigger than 20
                        if mode(s) < 0
                            Fx = - F_rep(abs(distance(j, s))); %repultion negative
                        else
                            Fx = F_att(abs(distance(j, s))); %attraction positive
                        end
                        total_force = total_force + Fx;
                    end
                end
            end
            % if total_force < 0
            %     flag = 1;
            %     break;
            % end 
            F_x(o) = total_force;
            o = o + 1;
            mover = mover + dx;
        end
        if flag == 0
        % Plot
        figure;
        plot(x_values, F_x);
        grid on;
        box on;
        xlabel('X-position');
        ylabel('Force applied to the system');
        title(['Force vs displacement with pitch ', num2str(k)]);
        end 
    end
end