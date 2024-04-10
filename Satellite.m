function [dY] = Satellite(t,Y,Mjd_UTC0,Thrust_LVLH,Xi)
%% Instantiate global variables
global AUXPAR EOPDATA CONST;

%% Demux
%Y = [x,y,z,x_dot,y_dot,z_dot,omega1,omega2,omega3,O_vec,Nu1;Nu2;Nu3]
PQR = Y(7:9);
O_vec = Y(10:18);
Nu = Y(19:21);

%% Orientation Matrices
O_EB = reshape(O_vec,3,3);

%% Compute Modified Julian time
Mjd_UTC = Mjd_UTC0+t/86400; % Coordinated Universal Time in julian time

% IERS dataset 
[x_pole,y_pole,UT1_UTC,LOD,~,~,~,~,TAI_UTC] = IERS(EOPDATA,Mjd_UTC,'l');

% Difference between UTC and TT
[~,~,~,TT_UTC,~] = TimeDiff(UT1_UTC,TAI_UTC); 

% Universal Time in julian time
Mjd_UT1 = Mjd_UTC+UT1_UTC/86400; 

% Terrestrial Time in julian time
Mjd_TT = Mjd_UTC+TT_UTC/86400; 

%% ICRS to ITRS transformation matrix and its derivative
P      = PrecMatrix(CONST.MJD_J2000,Mjd_TT);     % IAU 1976 Precession
N      = NutMatrix(Mjd_TT);                      % IAU 1980 Nutation
Theta  = GHAMatrix(Mjd_UT1,Mjd_TT);              % Earth rotation
Pi     = PoleMatrix(x_pole,y_pole);              % Polar motion
T = N*P;                                         % Transformation ECI to Truth-of-Date Inertia frame
E = Pi*Theta*T;                                  % Transformation ECI to ECEF

%% Gravitational perturbation from sun and moon
% Compute distance of sun and moon from earth with JPL DE440
MJD_TDB = Mjday_TDB(Mjd_TT);
[~,~,~,~,~,~,~,~,~,R_Moon,R_Sun,~] = JPL_Eph_DE440(MJD_TDB);

% Relative distance of sun and moon from satellite
R_SfrSC = R_Sun-Y(1:3);
R_MfrSC = R_Moon-Y(1:3);

% Compute gravitational perturbations from moon and sun
a_sun = CONST.GM_Sun*(R_SfrSC / norm(R_SfrSC)^3 - R_Sun/ norm(R_Sun)^3);
a_moon = CONST.GM_Moon*(R_MfrSC / norm(R_MfrSC)^3 - R_Moon/ norm(R_Moon)^3);

%% Earth's gravity model
if (AUXPAR.SolidEarthTides == 1 || AUXPAR.OceanTides == 1)
    %r_SC = norm(Y(1:3));
    %nominal_gravity = -CONST.GM_Earth*Y(1:3)/r_SC^3;
    %a_earth = nominal_gravity;
    a_earth = AccelHarmonic_AnelasticEarth(Mjd_UT1, Mjd_TT, R_Sun, R_Moon, Y(1:3),...
                    E, x_pole, y_pole);

else
    a_earth = AccelHarmonic(Y(1:3), E, AUXPAR.n, AUXPAR.m);

end

%% Thruster 
%Transform thrust from satellite body (LVLH) to Hill axis
Thrust_hill = [CONST.Thrust_map CONST.Thrust_map]*Thrust_LVLH;

% Transform to ECI axis
F_thrust = O_EB*Thrust_hill;

%% Solar radiation pressure
if (AUXPAR.SRAD == 1)
    % Cylindrical shadow model, no radiation if satellite behind
    % earth 
    e_Sun = R_Sun / norm(R_Sun);        % Sun direction unit vector
    s     = dot ( Y(1:3), R_Sun );      % Projection of s/c position 
    if ( s>0 || norm(Y(1:3)-s*e_Sun) > CONST.R_Earth )
        sun_radiate = 1;
    else
        sun_radiate = 0;
    end
    a_srad = sun_radiate*CONST.Csrp*CONST.area_solar*(1+CONST.Crefl)/2/CONST.Mass_sc...
                *R_SfrSC/norm(R_SfrSC);
else
    a_srad = zeros(3,1);

end

%% Atmospheric drag
if (AUXPAR.DRAG == 1)
	Omega = CONST.omega_Earth-0.843994809*1e-9*LOD; % IERS [rad/s]
    dens = NRLMSISE_00(Mjd_UTC,E*Y(1:3),UT1_UTC,TT_UTC);
    a_drag = AccelDrag(dens,Y(1:3),Y(4:6),T,CONST.area_drag,CONST.Mass_sc,CONST.Cd,Omega);

else
    a_drag = zeros(3,1);

end

%% Acceleration in F_E (ECI) frame
accel = a_earth + a_sun*0 + a_moon*0 + F_thrust/CONST.Mass_sc + a_srad + a_drag;
% [a_earth a_sun a_moon a_moon a_srad a_drag]  
dPos = [Y(4:6);accel];

%% Thruster torque 
tau_thrust = [CONST.Thrust_loc -CONST.Thrust_loc]*Thrust_LVLH; 

%% Reaction wheel dynamics
dNu = Xi;

%% Angular acceleration in F_L frame
dOmega = inv(CONST.J_sc)\(cross(CONST.J_sc*PQR + CONST.J_alp*Nu, PQR) - CONST.J_alp*Xi + tau_thrust);

%% Poisson Equation
crPQR = [0 -PQR(3) PQR(2);
    PQR(3) 0 -PQR(1);
    -PQR(2) PQR(1) 0];
dO_EB = -crPQR*O_EB;
dO_vec = dO_EB(:);


%% Output
dY = [dPos; dOmega; dO_vec; dNu];

end