%--------------------------------------------------------------------------
%
% MeanObliquity: Computes the mean obliquity of the ecliptic
%
% Input:
%  Mjd_TT    Modified Julian Date (Terrestrial Time)
% 
% Output:
%  MOblq     Mean obliquity of the ecliptic
%
% Last modified:   2018/01/27   Meysam Mahooti
% 
%--------------------------------------------------------------------------
function MOblq = MeanObliquity(Mjd_TT)

global CONST

T = (Mjd_TT-CONST.MJD_J2000)/36525;

MOblq = CONST.Rad*(23.43929111-(46.8150+(0.00059-0.001813*T)*T)*T/3600);

