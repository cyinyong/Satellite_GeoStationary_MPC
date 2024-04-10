%--------------------------------------------------------------------------
%
% ECEF2ECI: Transforms Earth Centered Earth Fixed (ECEF) coordinates to 
%           Earth Centered Inertial (ECI) coordinates
%
% Last modified:   2022/06/16   Meysam Mahooti
%
%--------------------------------------------------------------------------
function Y = ECEF2ECI(MJD_UTC, Y0)

global EOPDATA CONST

[x_pole,y_pole,UT1_UTC,LOD,~,~,~,~,TAI_UTC] = IERS(EOPDATA,MJD_UTC,'l');
[~,~,~,TT_UTC,~] = TimeDiff(UT1_UTC,TAI_UTC);

Y0 = reshape(Y0,[1,6]);
MJD_UT1 = MJD_UTC + UT1_UTC/86400;
MJD_TT  = MJD_UTC + TT_UTC/86400; 

% ICRS to ITRS transformation matrix and its derivative
P      = PrecMatrix(CONST.MJD_J2000,MJD_TT);     % IAU 1976 Precession
N      = NutMatrix(MJD_TT);                      % IAU 1980 Nutation
Theta  = GHAMatrix(MJD_UT1,MJD_TT);              % Earth rotation
Pi     = PoleMatrix(x_pole,y_pole);              % Polar motion

S = zeros(3);
S(1,2) = 1;
S(2,1) = -1;
% Omega = 7292115.8553e-11+4.3e-15*( (MJD_UTC-const.MJD_J2000)/36525 ); % [rad/s]
Omega  = CONST.omega_Earth-0.843994809*1e-9*LOD; % [rad/s] IERS
dTheta = Omega*S*Theta;                          % Derivative of Earth rotation matrix [1/s]
U      = Pi*Theta*N*P;                           % ICRS to ITRS transformation
dU     = Pi*dTheta*N*P;                          % Derivative [1/s]

% Transformation from WGS to ICRS
r = U'*Y0(1:3)';
v = U'*Y0(4:6)' + dU'*Y0(1:3)';
Y = [r;v];

