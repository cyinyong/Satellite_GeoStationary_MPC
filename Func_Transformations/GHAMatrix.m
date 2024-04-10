%--------------------------------------------------------------------------
%
% GHAMatrix: Transformation from true equator and equinox to Earth equator
%            and Greenwich meridian system 
%
% Inputs:
%   Mjd_UT1   Modified Julian Date UT1
%   Mjd_TT    Modified Julian Date TT
% 
% Output:
%   GHAmat    Greenwich Hour Angle matrix
%
% Last modified:   2022/06/16   Meysam Mahooti
%
%--------------------------------------------------------------------------
function GHAmat = GHAMatrix(Mjd_UT1,Mjd_TT)

GHAmat = R_z(GAST(Mjd_UT1,Mjd_TT));

