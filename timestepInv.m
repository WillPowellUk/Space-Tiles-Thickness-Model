function  [bestMethod, maxdt] = timestepInv(tmax, ntMin, ntMax, ntIncr, thick, nx, tol)

% [bestMethod, maxTimestep]
% Function for determining the appropriate timestep
% to maintain stability and accuracy, using the numerical
% approximation methods - forward, backward, Du-fort Frankel and 
% Crank-Nicolson

% W Powell  06/04/21
% Modified from code written by Dr Nigel Johnston

% Required input arguments:
% tmax   - maximum time
% ntMin, ntMax, ntIncr - minimum, maximum and incrementation of timesteps
%                        to test
% thick  - total thickness
% nx     - number of spatial steps in x axis
% tol    - error tolerance for accuracy

% Output arguments:
% bestMethod    - returns best method - method with maximum timestep 
%                 whilst maintaining accuracy and stability
% maxdt   - returns maximum timestep for best method, with optimal 
%           efficiency whilst remaing in the accepted error tolerance

% For example, to perform function with a max error of ±0.5 Kelvin:
% [bestMethod, maxdt] = timestepInv(4000, 41, 2001, 20, 0.05, 21, 0.5)


i = 0;
% iterate through values of nt for each method and extract final temp value
% at boundary
for nt = ntMin : ntIncr : ntMax
    i=i+1;
    dt(i) = tmax/(nt-1);
    disp (['nt = ' num2str(nt) ', dt = ' num2str(dt(i)) ' s']);
    [~, ~, u] = shuttle(tmax, nt, thick, nx, 'forward', false);
    uf(i) = u(end, 1);
    [~, ~, u] = shuttle(tmax, nt, thick, nx, 'backward', false);
    ub(i) = u(end, 1);
    [~, ~, u] = shuttle(tmax, nt, thick, nx, 'dufort-frankel', false);
    udf(i) = u(end, 1);
    [~, ~, u] = shuttle(tmax, nt, thick, nx, 'crank-nicolson', false);
    ucn(i) = u(end, 1);
end
% calculates average stable value of temperature for all method at smallest
% dt value
avgConvergence = (uf(end) + ub(end) + udf(end) + ucn(end))/4;
plot(dt, [uf; ub; udf; ucn])
hold on
yline(avgConvergence + tol, 'k--')
yline(avgConvergence - tol, 'k--')
hold off
ylim([avgConvergence - 3*tol, avgConvergence + 3*tol])
xlabel('Timestep, dt (s)')
ylabel('Inner tile temperature at 4000 seconds (K)')
legend ('Forward', 'Backward', 'Dufort-frankel', 'Crank-nicolson')

maxdt = 0;
i = 0;
methods = {'Forward', 'Backward', 'Dufort-frankel', 'Crank-nicolson'};
z = [uf; ub; udf; ucn];
% iterates through methods and finds method with largest timestep whilst
% still being in the tolerance range
for i = 1:length(methods)
    method = methods(i);
    intersectPos = [];
    intersectNeg = [];
    intersectPos = dt(find(z(i,:) > avgConvergence + tol, 1, 'last')+1);
    intersectNeg = dt(find(z(i,:) < avgConvergence - tol, 1, 'last')+1);
    if ~isempty(intersectPos) && intersectPos > maxdt
        maxdt = intersectPos;
        bestMethod = method;
    elseif ~isempty(intersectNeg) && intersectNeg > maxdt
        maxdt = intersectNeg;
        bestMethod = method;
    end
end

end

