%% Set global constants and read Models
global EOPDATA PC Cnm Snm SWDATA CONST

%% Constants
% Mathematical constants
CONST.pi2       = 2*pi;                          % 2pi
CONST.Rad       = pi/180;                        % Radians per degree
CONST.Deg       = 180/pi;                        % Degrees per radian
CONST.Secs      = 86400;                         % Seconds per day
CONST.Arcs      = 3600*180/pi;                   % Arcseconds per radian
CONST.TURNAS    = 1296000.0;                     % Arcseconds in a full circle
CONST.DAS2R     = 4.848136811095359935899141e-6; % Arcseconds to radians

% General
CONST.MJD_J2000 = 51544.5;             % Modified Julian Date of J2000
CONST.T_B1950   = -0.500002108;        % Epoch B1950
CONST.c_light   = 299792457.999999984; % Speed of light  [m/s]; DE440
CONST.AU        = 149597870699.999988; % Astronomical unit [m]; DE440

% Physical parameters of the Earth, Sun and Moon

% Equatorial radius and flattening
CONST.R_Earth   = 6378.137e3;      % Earth's radius [m]; WGS-84
CONST.f_Earth   = 1/298.257223563; % Flattening; WGS-84
CONST.R_Sun     = 696000.0e3;      % Sun's radius [m]; DE440
CONST.R_Moon    = 1738.0e3;        % Moon's radius [m]; DE440

% Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
CONST.omega_Earth = 15.04106717866910/3600*CONST.Rad; % [rad/s]; WGS-84

% Gravitational coefficients
CONST.GM_Earth   = 398600.4415e9;                      % [m^3/s^2]; GGM03C & GGM03S
CONST.GM_Sun     = 132712440041.279419e9;              % [m^3/s^2]; DE440
CONST.GM_Moon    = CONST.GM_Earth/81.3005682214972154; % [m^3/s^2]; DE440
CONST.GM_Mercury = 22031.868551e9; 		               % [m^3/s^2]; DE440
CONST.GM_Venus   = 324858.592000e9;                    % [m^3/s^2]; DE440
CONST.GM_Mars    = 42828.375816e9;	      			   % [m^3/s^2]; DE440
CONST.GM_Jupiter = 126712764.100000e9;    			   % [m^3/s^2]; DE440
CONST.GM_Saturn  = 37940584.841800e9;     			   % [m^3/s^2]; DE440
CONST.GM_Uranus  = 5794556.400000e9;      			   % [m^3/s^2]; DE440
CONST.GM_Neptune = 6836527.100580e9;      			   % [m^3/s^2]; DE440
CONST.GM_Pluto   = 975.500000e9;	      			   % [m^3/s^2]; DE440

% Solar radiation pressure at 1 AU 
%CONST.P_Sol = 1367/CONST.c_light; % [N/m^2] (1367 W/m^2); IERS 96

% Satellite Parameters
CONST.Mass_sc = 4000;                             % mass of satellite
CONST.J_sc = diag([1.7e4;2.7e4;2.7e4]);           % MOI of satellite bus
CONST.J_alp = eye(3)*0.8;                         % MOI of reaction wheels
CONST.Thrust_loc = [0 2.5 0; 0 0 2.5; 3.75 0 0];  % Dual Thruster location from CG 
CONST.Thrust_map = [0 0 1; 1 0 0; 0 1 0];         % Dual Thruster mapping

% Solar radiation
CONST.Csrp = 9.1e-6;      % Solar radiation pressure coefficient
CONST.Crefl = 0.6;        % Reflection coefficient
CONST.area_solar = 200;   % Solar radiation area

% Wind Drag
CONST.area_drag = 7.5*5;  % Wind drag area                  
CONST.Cd = 2.7;           % Wind drag coefficient

%% Load DE440Coeff 
load DE440Coeff.mat
PC = DE440Coeff;

%% Load Earth gravity field coefficients
Cnm = zeros(361,361);
Snm = zeros(361,361);
fid = fopen('GGM03C.txt','r');
for n=0:360
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end
fclose(fid);

%% Load Earth orientation parameters into EOPDATA
fid = fopen('EOP-All.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
while ~feof(fid)
    tline = fgetl(fid);
    k = strfind(tline,'NUM_OBSERVED_POINTS');
    if (k == 1)
        numrecsobs = str2num(tline(21:end));
        tline = fgetl(fid);
        for i=1:numrecsobs
            EOPDATA(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        for i=1:4
            tline = fgetl(fid);
        end
        numrecspred = str2num(tline(22:end));
        tline = fgetl(fid);
        for i=numrecsobs+1:numrecsobs+numrecspred
            EOPDATA(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        break
    end
end
fclose(fid);

%% Load space weather data into SWDATA
fid = fopen('SW-All.txt','r');
%  ---------------------------------------------------------------------------------------------------------------------------------
% |                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
% | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
%  ---------------------------------------------------------------------------------------------------------------------------------
while ~feof(fid)
    tline = fgetl(fid);
    k = strfind(tline,'NUM_OBSERVED_POINTS');
    if (k == 1) 
        numrecsobs = str2num(tline(21:end));
        tline = fgetl(fid);
        for i=1:numrecsobs
            SWDATA(:,i) = fscanf(fid,'%i %d %d %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %f %i %i %f %i %f %f %f %f %f',[33 1]);
        end
        % remove the row of the Q parameter
        SWDATA = [SWDATA(1:27,:);SWDATA(29:33,:)];
        for i=1:4
            tline = fgetl(fid);
        end
        numrecspred = str2num(tline(28:end));
        tline = fgetl(fid);
        %  -------------------------------------------------------------------------------------------------------------------------------
        % |                                                                                             Adj   Adj   Adj   Obs   Obs   Obs 
        % | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Ctr81 Lst81 F10.7 Ctr81 Lst81
        %  -------------------------------------------------------------------------------------------------------------------------------
        for i=numrecsobs+1:numrecsobs+numrecspred
            SWDATA(:,i) = fscanf(fid,'%i %d %d %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %f %i %i %f %f %f %f %f %f',[32 1]);
        end
        break
    end
end
fclose(fid);