function [x, t, u, uEnd] = shuttle2d(tmax, nt, xmax, nx, ymax, ny, method, varargin)

% Function for modelling temperature in a space shuttle tile
% W Powell  06/04/2021
% Co-Author: F Berteau

% Required input arguments:
% tmax   - maximum time
% nt     - number of timesteps
% xmax   - total thickness of tile
% nx     - number of spatial steps in x axis
% ymax   - width of tile
% ny     - number of spatial steps in y axis
% method - solution method ('forward', 'backward' etc)

% Optional input arguments:
% doPlot - true to plot graph; false to suppress graph.
% tempUnitK - true if Kelvin, false if Fahrenheit
% tileMat - tile material
% fileName - graph image for data extraction
% axes - axes handle to plot to GUI; false to plot to window

% Return arguments:
% x      - distance vector
% t      - time vector
% u      - temperature matrix
%
% For example, to perform a 2d simulation for 1001 time steps:
%   [x, t, u] = shuttle2d(4000, 1001, 0.05, 21, 0.2, 81, 'forward');


% optional parameters
switch nargin
    case {1,2,3,4,5, 6, 7} 
        doPlot = true;
        tempUnitK = true;
        tileMat = 'li900';
        fileName = '597';
        axes = false;
    case 8
        doPlot = varargin{1};
        tempUnitK = true;
        tileMat = 'li900';
        fileName = '597';
        axes = false;
    case 9
        [doPlot, tempUnitK] = varargin{1:2};
        tileMat = 'li900';
        fileName = '597';
        axes = false;
    case 10
        [doPlot, tempUnitK, tileMat] = varargin{1:3};
        fileName = '597';
        axes = false;
    case 11
        [doPlot, tempUnitK, tileMat, fileName] = varargin{1:4};
        axes = false;
    case 12
        [doPlot, tempUnitK, tileMat, fileName, axes] = varargin{1:5};
end

% Set tile properties
switch tileMat
    case 'li900'
        thermCon = 0.0577; % W/(m K)
        density  = 144;   % 9 lb/ft^3
        specHeat = 1261;  % ~0.3 Btu/lb/F at 500F
    case 'li2200'
        thermCon = 0.0577; % W/(m K)
        density  = 352.5;   % 22 lb/ft^3
        specHeat = 1261;  % ~0.3 Btu/lb/F at 500F
    case 'frcp12'
        thermCon = 0.0577; % W/(m K)
        density  = 192.2;   % 12 lb/ft^3
        specHeat = 1261;  % ~0.3 Btu/lb/F at 500F
end

% loads data from graph image using image extraction function
[tempK, tempF, time] = imgExtraction(fileName);

% Initialise everything.
dt = tmax / (nt-1); % sets time step size
t = (0:nt-1) * dt; % time vector
dx = xmax / (nx-1); % sets thickness step size
x = (0:nx-1) * dx; % thickness of tile as vector
dy = ymax / (ny-1); % tile width step size
y = (0:ny-1) * dy; % tile width as vector
u = zeros(ny, nx, nt); % pre allocates space for U (state vector for temperature)
alpha = thermCon/(density * specHeat); % thermal diffusivity
p = (alpha * dt)/(dx^2); % non dimensional timestep

% Use interpolation to get outside temperature at times t 
% and store it as right-hand boundary R.
if tempUnitK
    R = interp1(time, tempK, t, 'linear', tempK(end));
else
    R = interp1(time, tempF, t, 'linear', tempF(end));
end

% set initial conditions equal to boundary temperature at t=0.
u(:, :, 1) = R(1);

% plots a animated graph. h is a 'handle' of the graph
if doPlot
    if axes == false
        figure(3);
        h = surf(x, y, u(:,:,1)); % surf plot
        pbaspect([1 4 1]) % sets the aspect ratio
        %displays the current time of animation
        h2=text(-0.05, 0.18, 280, ' t = 0 s   ', 'BackgroundColor', [1,1,1]);
        shading interp % shades for clearer graphics
        % set y axis range dependent on Kelvin or Fahrenheit
        if tempUnitK
            tempUnit = 'K';
            zlim([280 1300])
        else
            tempUnit = 'F';
            zlim([44 1880])
        end
        % colour table to translate colour to a temperature value
        caxis ([0 1200]);
        colorbar
        % set axes labels
        xlabel('\itx\rm (m)')
        ylabel('\ity\rm (m)')
        zlabel(['\itu\rm (' tempUnit ')'])
    else
        h = surf(axes, x, y, u(:,:,1)); 
        set(h, 'EdgeColor', 'interp', 'FaceColor', 'interp'); % shades for clearer graphics
        h2=text(axes, -0.03, 0.13, 280, ' t = 0 s   ', 'BackgroundColor', [1,1,1]);
        % colour table to translate colour to a temperature value
        caxis (axes,[0 1200]); 
        colorbar(axes) 
    end
end

% Main timestepping loop.
for n = 1:nt-1

    % Select method.
    switch method
        % forward method with Neuman Boundary Condition
        case 'forward'
            % Boundary condition from outside temperature
            u(:,nx,n+1) = R(n+1);
            
            % Neuman Boundaries
            i = 1:nx-1;
            im = [2 1:nx-2];
            ip = 2:nx;

            j = 2:ny;
            jm = 1:ny-1;
            jp = [3:ny ny-1];

            % calculates U state vector with respect to thickness and width of the tile
            u(j, i, n+1) = (1 - 4 * p) * u(j, i, n) + ...
                p * (u(j, im, n) + u(j, ip, n) + u(jm, i, n) + u(jp, i, n));
            u(1, i, n+1) = (1 - 4 * p) * u(1, i, n) + ...
                p * (u(1, im, n) + u(1, ip, n) + 2*u(2, i, n));
        % backward method without Neuman Boundary
        case 'backward'
            maxiterations = 100; % timestep
            tolerance = 1.e-4; % for stability
            % Boundary condition from outside temperature
            u(:,nx,:) = R(n+1);
            % calculate internal values iteratively using Gauss-Seidel
            % Starting values are equal to old values
            u(:, :, n+1) = u(:, :, n);

            for iteration = 1:maxiterations
                change = 0;
                for i=2:nx-1
                    for j=2:ny-1
                        uold = u(j, i, n+1);
                        u(j, i, n+1) = ((u(j, i, n) + p * (u(j-1, i, n+1) + u(j+1, i, n+1)...
                            + u(j, i-1, n+1) + u(j, i+1, n+1)))/(1+4*p));
                        change = change + abs(u(j, i, n+1) - uold);
                    end
                end
                % exit backward method calculations if unstable
                if change < tolerance
                    break
                end
            end
            %disp(['Time = ' num2str(t(n)) ' s: Iterations = ' num2str(iteration)]);
        otherwise
            error (['Undefined method: ' method])
            return
    end
    % update graph with new values
    set(h,'ZData', u(:, :, n+1));
            
    % display current time
    txt = sprintf(' t = %4.1f s ', t(n+1));
    set(h2, 'String', txt)
    % Refresh screen
    drawnow
end
uEnd = u(end);
end


    