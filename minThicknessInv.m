function [thickness] = minThicknessInv(maxTemp, tmax, nt, nx, method, tempUnitK, doPlot)

% This function determines the minimum tile thickness required to prevent
% an overshoot of an inputted maximum temperature at the inside of the tile.

% W Powell  06/04/21

% Required input arguments:
% maxTemp - maximum temperature desired at inside of tile
% tmax   - maximum time
% nt     - number of timesteps
% xmax   - total thickness
% nx     - number of spatial steps in x axis
% method - solution method ('forward', 'backward' etc)
% tempUnitK - true if Kelvin, false if Fahrenheit
% doPlot - true to plot graph; false to suppress graph.

% Output arguments:
% thickness    - minimum tile thickness required to withold from maxTemp

% For example, to perform function for a maximim temperatre of 450K:
% [thickness] = minThicknessInv(450, 4000, 501, 21, 'crank-nicolson', true, true)


% initial guesses for minimum thickness
guess(1,1) = 0.05;
guess(1,2) = 0.1;

% Desired temperature lower than maximum temperature(K or F)
tempPrecision = 1;

% first guess for shooting method
[~, ~, u] = shuttle(tmax, nt, guess(1,1), nx, method, false);
guess(2,1) = max(u(:,1)); % finds maximum temperature

% second guess
[~, ~, u] = shuttle(tmax, nt, guess(2,1), nx, method, false);

% sets n to third guess
n = 3;

% Shooting Method will continue to iterate until the final guess has
% acheived the desired precision for maximum temperature of the tile
while guess(2,end) > maxTemp || guess(2,end) < (maxTemp - tempPrecision)
    
    % work out the linear gradient, m, of the last two guesses
    m = (guess(2, n-2) - guess(2, n-1)) / (guess(1,n-2) - guess(1,n-1));
    
    % find the y intercept, c, using m
    c = guess(2, n-1) - m * guess(1, n-1);
    
    % An intelligent new guess of tile thickness is added to the nth value 
    % of tile thickness guesses. The next guess is computed using 
    % the Secant root finding method
    guess(1, n) = (maxTemp - c) / m;
    
    % finds maximum temperature in either K or F at that tile thickness
    if tempUnitK
        [~, ~, u] = shuttle(tmax, nt, guess(1,n), nx, method, false);
        guess(2, n) = max(u(:,1));
    else
        [~, ~, u] = shuttle(tmax, nt, guess(1,n), nx, method, false, false);
        guess(2, n) = max(u(:,1));
    end
    % prepares for next guess 
    n = n + 1;
end

% final guess is used for minimum thickness
thickness = guess(1, end);

% plots 3d graph of tile temperature with respect to tile thickness, with time domain
if doPlot
    if tempUnitK
        shuttle(tmax, nt, guess(1,end), nx, method, true);
    else
        shuttle(tmax, nt, guess(1,end), nx, method, true, false);
    end
end

end