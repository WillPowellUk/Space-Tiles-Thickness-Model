function [x, t, u] = shuttle(tmax, nt, xmax, nx, method, varargin)

% Function for modelling temperature in a space shuttle tile

% W Powell  06/04/21
% Modified from code written by Dr Nigel Johnston

% Required input arguments:
% tmax   - maximum time
% nt     - number of timesteps
% xmax   - total thickness
% nx     - number of spatial steps in x axis
% method - solution method ('forward', 'backward' etc)

% Optional input arguments:
% doPlot - true to plot graph; false to suppress graph.
% tempUnitK - true if Kelvin, false if Fahrenheit
% tileMat - tile material
% fileName - graph image for data extraction
% GUI - plot to GUI; false to plot to window
% axes1, axes2, axes3 - handle for axes (only for GUI)

% Return arguments:
% x      - distance vector
% t      - time vector
% u      - temperature matrix
%
% For example, to perform a  simulation with 501 time steps
%   [x, t, u] = shuttle(4000, 501, 0.05, 21, 'forward', true);



% optional parameters
switch nargin
    case {1,2,3,4,5} 
        doPlot = true;
        tempUnitK = true;
        tileMat = 'li900';
        fileName = '597';
        GUI = false;
    case 6
        doPlot = varargin{1};
        tempUnitK = true;
        tileMat = 'li900';
        fileName = '597';
        GUI = false;
    case 7
        [doPlot, tempUnitK] = varargin{1:2};
        tileMat = 'li900';
        fileName = '597';
        GUI = false;
    case 8
        [doPlot, tempUnitK, tileMat] = varargin{1:3};
        fileName = '597';
        GUI = false;
    case 9
        [doPlot, tempUnitK, tileMat, fileName] = varargin{1:4};
        GUI = false;
    case 10
        [doPlot, tempUnitK, tileMat, fileName, GUI] = varargin{1:5};
    case 13
        [doPlot, tempUnitK, tileMat, fileName, GUI, axes1, axes2, axes3 ] = varargin{1:8};
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
        specHeat = 1261;
        % ~0.3 Btu/lb/F at 500F
end

% loads data from graph image using image extraction function
[tempK, tempF, time] = imgExtraction(fileName);

% Initialise everything.
dt = tmax / (nt-1); % sets time step size
t = (0:nt-1) * dt; % time vector
dx = xmax / (nx-1); % sets thickness step size
x = (0:nx-1) * dx; % thickness of tile as vector
u = zeros(nt, nx); % pre allocates space for U (state vector for temperature)
alpha = thermCon/(density * specHeat); % thermal diffusivity
p = (alpha * dt)/(dx^2); % non dimensional timestep

% Use interpolation to get outside temperature at times t 
% and store it as right-hand boundary R.
if tempUnitK
    tempUnit = 'K';
    R = interp1(time, tempK, t, 'linear', tempK(end));
else
    tempUnit = 'F';
    R = interp1(time, tempF, t, 'linear', tempF(end));
end
% set initial conditions equal to boundary temperature at t=0.
u(1, :) = R(1);

% Main timestepping loop.
for n = 1:nt-1

    % Select method.
    switch method
        case 'forward'
            % neuman boundary conditions
            i = 1:nx-1;
            im = [2 1:nx-2];
            ip = 2:nx;
            % set boundary conditions
            u(n+1,nx) = R(n+1);
            u(n+1,i) = (1 - 2*p) * u(n,i) + p * (u(n,im) + u(n,ip));
        case 'dufort-frankel'
            % neuman boundary conditions
            i=1:nx-1;
            im = [2 1:nx-2];
            ip = 2:nx;
            % set boundary conditions
            u(n+1,1) = R(1);
            u(n+1,nx) = R(n+1);
            % set index for 'old' point
            if n == 1
                nminus1 = 1; % at first timestep, old point doesn't exist as n-1 = 0.
                             % Use value at timestep 1 instead.
            else
                nminus1 = n-1; % after first timestep, proceed normally.
            end
            % calculate internal values using Leapfrog method
            u(n+1,i)= (u(nminus1,i) * (1-2*p) + 2 * p * (u(n,im) + u(n,ip)))/(1+2*p);
        case 'backward'
            ivec = 2:nx-1; % set up index vector
            % set boundary conditions
            l = R(1);
            r = R(n+1);
            ivec = 2:nx-1; % set up index vector
            %calculate internal values using backward differencing
            b(1)    = 1 + 2*p;
            c(1)    = -2*p;
            d(1)    = u(n,1);
            a(ivec) = -p;
            b(ivec) = 1 + 2*p;
            c(ivec) = -p;
            d(ivec) = u(n,ivec);
            a(nx)   = 0;
            b(nx)   = 1;
            d(nx)   = r;
            u(n+1,:) = triDiag(a,b,c,d); % calls tri diagonal matrix function
        case 'crank-nicolson'
            ivec = 2:nx-1; % set up index vector
            % set boundary conditions
            l = R(1);
            r = R(n+1);
            ivec = 2:nx-1; % set up index vector
            % calculate internal values using crank nicolson
            b(1)    = 1 + 2*p;
            c(1)    = -2*p;
            d(1)    = u(n,1);
            a(ivec) = -p/2;
            b(ivec) = 1 + p;
            c(ivec) = -p/2;
            d(ivec) = p/2*u(n,ivec-1) + (1-p)*u(n,ivec) + p/2*u(n,ivec+1);
            a(nx)   = 0;
            b(nx)   = 1;
            d(nx)   = r;
            u(n+1,:) = triDiag(a,b,c,d); % calls tri diagonal matrix function
        otherwise
            error (['Undefined method: ' method])
            return
    end
end

InsideTemp = u(:, 1);

% plots 3d graph of tile temperature with respect to tile thickness, with time domain
if doPlot
    if GUI == false % if axes location is undefined plot normally with window
        % plots inside and outside temperature of tile in the time domain
        if tempUnitK
            figure(1);
            plot(time, tempK);
            xlabel('Time (s)');
            ylabel('Temperature (K)');
            title('Outside Temperature of Tile');
            figure(2);
            plot(t, InsideTemp');
            xlabel('Time (s)');
            ylabel('Temperature (K)');
            title('Inside Temperature of Tile');
            xlim([0, tmax]);
            
        else
            figure(1);
            plot(time, tempF);
            xlabel('Time (s)');
            ylabel('Temperature (F)');
            title('Outside Temperature of Tile');
            figure(2);
            plot(t, InsideTemp');
            xlabel('Time (s)');
            ylabel('Temperature (F)');
            title('Inside Temperature of Tile');
            xlim([0, tmax]);
        end
        figure(3); 
        waterfall(x,t,u); % waterfall 3d graph plot
    else % if using GUI, plot in GUI not window
        % plots inside and outside temperature of tile in the time domain
        if tempUnitK
            plot(axes1, time, tempK);
            plot(axes2, t, InsideTemp');
        else
            plot(axes1, time, tempF);
            plot(axes2, t, InsideTemp');
        end
        waterfall(axes3,x,t,u); % waterfall 3d graph plot
    end
    view(150,30) % viewing angle of plot
    xlabel('\itx\rm (m)')
    ylabel('\itt\rm (s)')
    zlabel(['\itu\rm (' tempUnit ')'])

end
end



    