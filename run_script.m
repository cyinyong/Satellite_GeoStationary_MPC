clear all; close all; clc;

% Path settings
path(path,".\Func_Earth_Models");
path(path,".\Func_Mjd_Time");
path(path,".\Func_Transformations");

% Plot settings
set(0,'defaultAxesFontName', 'times',...
    'defaultTextFontName','times',...
    'DefaultAxesFontSize', 14, ...
    'DefaultTextFontSize', 14, ...
    'DefaultLineLineWidth',1.5,...
    'DefaultLineMarkerSize',6,...
    'DefaultTextInterpreter','latex')

% Initialise parameters
set_model_params;
global AUXPAR;

%% Simulation

% Geo stationary orbit
Total_Steps = 1000;                      % Total simulation step
Ts = 6*60*60;                           % Sampling time [s]

LLA0 = [0 40 35786e3];                  % Initial Lat Lon Alt [deg deg m]
Vel0 = [0;0;0];                         % Initial ECEF velocity (m/s) (assume geostationary orbit) 
Eul0 = [0;0;0];                         % Initial Euler angle (rad)
PQR0 = [0;0;0];                         % Initial body rates (rad/s)

% Low earth orbit
% Total_Steps = 1000;                      % Total simulation step
% Ts = 60; %2*60*60;                           % Sampling time [s]
% 
% LLA0 = [-4.07697364592954;
%     1.74513198492441;	
%     788049.103010912]';                  % Initial Lat Lon Alt [deg deg m]
% Vel0 = [562.650611; 
%     -1616.516697; 
%     7358.157263];                       % Initial ECEF velocity (m/s)
% Eul0 = [0;0;0];                         % Initial Euler angle (rad)
% PQR0 = [0;0;0];                         % Initial body rates (rad/s)

year = 2002;
month = 4;
day = 24;
hour = 21;
minute = 55;
sec = 28;
Mjd_UTC0 = Mjday(year, month, day, hour, minute, sec);

Omat0 = Eul2Omat(Eul0);
O_vec0 = zeros(9,1);
Pos_ECEF0 = [lla2ecef(LLA0,'WGS84')'; Vel0]; 
Pos_ECI0 = ECEF2ECI(Mjd_UTC0, Pos_ECEF0);

%% Set AUXPAR
% Solar radiation pressure
AUXPAR.SRAD = 1;

% Atmospheric drag
AUXPAR.DRAG = 1;

% Gravity models
AUXPAR.SolidEarthTides = 1; 
AUXPAR.OceanTides = 1;
AUXPAR.n       = 70; 
AUXPAR.m       = 70;

%% Simulation
% Set output matrix
y1 = zeros(21,Total_Steps);
Thrust_in = zeros(6,1);
Xi_in = zeros(3,1);
Pos_ECEF = zeros(6,Total_Steps);
LLA = zeros(3,Total_Steps);
Eul = zeros(3,Total_Steps);

% Initial state vector
y0_dyn = [Pos_ECI0; PQR0; O_vec0; Xi_in];

for ii=1:Total_Steps
    tic

    % Julian time
    Mjd_UTC = Mjd_UTC0+ii*Ts/86400;

    % Do control
    %Thrust_in = Thrust_in;
    %Xi_in = Xi_in

    % Propagate dynamic states
    if ii == 1
        [~,yy] = ode45(@(t,y)Satellite(t,y,Mjd_UTC,Thrust_in,Xi_in),[ii ii+1]*Ts, y0_dyn);
    else
        [~,yy] = ode45(@(t,y)Satellite(t,y,Mjd_UTC,Thrust_in,Xi_in),[ii ii+1]*Ts, y1(:,ii-1));
    end
    
    % Store states
    y1(:,ii) = yy(end,:);
    
    % Output
    Pos_ECEF(:,ii) = ECI2ECEF(Mjd_UTC, y1(1:6,ii));
    LLA(:,ii) = ecef2lla(Pos_ECEF(1:3,ii)','WGS84')';
    Omat = reshape(y1(10:18,ii),[3,3]);
    Eul(:,ii) = Omat2Eul(Omat);
    disp(strcat('iter:',mat2str(ii)))
    
    toc
end

% Add initial condition 
y1 = [y0_dyn y1];
LLA = [LLA0' LLA];
Eul = [Eul0 Eul];

%% Plot
time = [0:1:Total_Steps] * Ts / 60/60;

h = figure(1);
h.Color = 'w';
h.Units = 'inches';
ha = tight_subplot(1,2,[.1 .1],[.12 .05],[.12 .05]);

axes(ha(1))
plot(LLA(2,:),LLA(1,:))
grid
ylabel('Latitude error $(deg)$')
xlabel('Longitude error $(deg)$')

axes(ha(2))
plot(time,LLA(3,:)/10^3)
grid
ylabel('Alt$(km)$')
xlabel('t$(hr)$')

h=figure(2);
h.Color = 'w';
plot3(y1(1,:)/10^3,y1(2,:)/10^3,y1(3,:)/10^3)
grid
xlabel('$X(km)$')
ylabel('$Y(km)$')
zlabel('$Z(km)$')

%%
figure(4)
geoshow('landareas.shp','FaceColor',[0.5 1 0.5]);
title('Satellite''s Ground Track')
hold on
plot(LLA(2,:),LLA(1,:),'.r','MarkerSize',15)
% animation
an = animatedline('Marker','*','MarkerSize',25);
for k = 1:Total_Steps+1
    addpoints(an,LLA(2,k),LLA(1,k));
    drawnow
    pause(0.01);
    clearpoints(an);
end

%%
uif = uifigure;
g = geoglobe(uif);
geoplot3(g,LLA(1,:),LLA(2,:),LLA(3,:),"r")


