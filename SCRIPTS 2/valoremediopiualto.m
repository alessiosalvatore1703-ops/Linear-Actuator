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

% Part that interpolates the values and finds an approximation of the force
F_att = fit(d, FxAttractioncut, 'cubicspline');
F_rep = fit(d, FxRepultioncut, 'cubicspline');

Values_attraction = F_att(test_values);
Values_repultion = F_rep(test_values);

F_x = zeros(1, length(x_values)); 
num_magn = 4;

% Initialize a matrix to store the force values and corresponding k, p, f combinations
force_values = [];

for k = 22 : 1 : 22  % Minimum pitch equal to coil pitch, max double
    for p = k+15 : 1 : 40
        for f = 125 : 1 : 137
            mover = [0, k , p, f];
            coils = pitch * (0 : 1 : num_coils);
            o = 1;
            flag = 0;
            total_forces = 0; % Initialize total force sum for averaging
            num_iterations = 0;

            % Simulate the interaction between the magnet and the coils and
            % store the force value
            for x = x_values
                distance = zeros(num_magn, num_coils); % Matrix storing the distance of each magnet to each coil
                mode = zeros(1, num_coils); % Array containing the mode of the coil
                for s = 1 : num_coils  
                    for j = 1 : num_magn
                         distance(j, s) = coils(s) - mover(j);
                    end
                end

                sum = 0;
                total_force = 0;

                % Define the mode of the coils for each iteration
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
                        if abs(distance(j, s)) < 20 % No information if greater than 20
                            if mode(s) < 0 
                                Fx = - F_rep(abs(distance(j, s))); % Repulsion negative
                            else 
                                Fx = F_att(abs(distance(j, s))); % Attraction positive
                            end
                            total_force = total_force + Fx;
                        end
                    end
                end

                if total_force < 0
                    flag = 1;
                    break;
                end

                F_x(o) = total_force;
                o = o + 1;
                mover = mover + dx;

                % Accumulate forces for averaging
                total_forces = total_forces + total_force;
                num_iterations = num_iterations + 1;
            end

            if flag == 0
                avg_force = total_forces / num_iterations;  % Calculate average force for this combination

                % Store the combination and its average force
                force_values = [force_values; k, p, f, avg_force];
            end
        end
    end
end

if size(force_values, 1) >= 10
    % Sort the combinations by the average force in descending order
    sorted_force_values = sortrows(force_values, 4, 'descend');

    % Get the top 10 combinations with the highest average force
    top_10_combinations = sorted_force_values(1:10, :);

    % Display the top 10 combinations and their average force
    disp('Top 10 Combinations with Highest Average Force:');
    disp('k    p    f    Avg Force');
    disp(top_10_combinations);
else
    disp('There are less than 10 valid combinations.');
    disp('Displaying all available combinations:');
    disp(force_values);
end