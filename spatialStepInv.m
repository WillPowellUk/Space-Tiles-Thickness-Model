function  [bestMethod, maxdx] = spatialStepInv(tmax, nxMin, nxMax, nxIncr, thick, nt, tol)

% Function for determining the appropriate Spatial Step
% to maintain stability and accuracy, using the numerical
% approximation methods - forward, backward, Du-fort Frankel and 
% Crank-Nicolson

% W Powell  06/04/21
% Modified from code written by Dr Nigel Johnston

% Required input arguments:
% tmax   - maximum time
% nxMin, nxMax, nxIncr - minimum, maximum and incrementation of timesteps
%                        to test
% thick  - total thickness
% nt     - number of timesteps, use stable value for all methods 
%          from timestepInv
% tol    - error tolerance for accuracy

% Output arguments:
% bestMethod    - returns best method - method with maximum timestep 
%                 whilst maintaining accuracy and stability
% maxdx   - returns maximum spatial step for best method, with optimal 
%           efficiency whilst remaing in the accepted error tolerance

% For example, to perform function with a max error of ±0.5 Kelvin:
%[bestMethod, maxdx] = spatialStepInv(4000, 5, 30, 1, 0.05, 1001, 0.5)


i = 0;
% iterate through values of nx for each method and extract final temp value
% at boundary
for nx = nxMin : nxIncr : nxMax
    i=i+1;
    dx(i) = thick / (nx-1);
    disp (['nx = ' num2str(nx) ', dx = ' num2str(dx(i)*1000) 'mm']);
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
% avgConvergence = (uf(end) + ub(end) + udf(end) + ucn(end))/4;
avgConvergence = 431.962935806185;
plot(dx*1000, [uf; ub; udf; ucn])
hold on
yline(avgConvergence + tol, 'k--')
yline(avgConvergence - tol, 'k--')
hold off
ylim([avgConvergence - 3*tol, avgConvergence + 3*tol])
xlabel('Spatialstep, dx (mm)')
ylabel('Inner tile temperature at 4000 seconds (K)')
legend ('Forward', 'Backward', 'Dufort-frankel', 'Crank-nicolson')

maxdx = 0;
i = 0;
methods = {'Forward', 'Backward', 'Dufort-frankel', 'Crank-nicolson'};
z = [uf; ub; udf; ucn];
% iterates through methods and finds method with largest timestep whilst
% still being in the tolerance range
% for i = 1:length(methods)
%     method = methods(i);
%     intersectPos = [];
%     intersectNeg = [];
%     intersectPos = dx(find(z(i,:) > avgConvergence + tol, 1, 'last')+1);
%     intersectNeg = dx(find(z(i,:) < avgConvergence - tol, 1, 'last')+1);
%     if ~isempty(intersectPos) && intersectPos > maxdx
%         maxdx = intersectPos;
%         bestMethod = method;
%     elseif ~isempty(intersectNeg) && intersectNeg > maxdx
%         maxdx = intersectNeg;
%         bestMethod = method;
%     end
% end

end

