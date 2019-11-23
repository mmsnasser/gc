clear
%
%
clc
tic
n      =  2^4;
t      = (0:2*pi/n:2*pi-2*pi/n).';
ratio  =  1;
%
% 
level =  4;
number_of_intervales = 2^level
c = 0; r = 1.5; % level of finite Cantor set
for k = 1:level
r = r/3; c = [c+2*r c-2*r];
end
Lc     =  sort(c.');
Lk     =  1/(3^(level-1))+zeros(size(Lc));
thetk  =  zeros(size(Lc));
%
%
%
[et , etp , cent , fet] = PreImageStrSlit (Lc , Lk , thetk , ratio , n , 1e-14 , 100 );
zet  =  et+fet;
%
m        =  length(Lc);
alphav   =  cent;
deltav   =  zeros(size(alphav)); deltav(m/4+1:3*m/4)=1;
mp       =   m; 
ell      =   0;
%
%
%
% 
%
[cap,u0] = capgc (et,etp,alphav,deltav,m,mp,ell,inf,0);
har_meas_0_ie_method = u0;
time_ie_method = toc
%
har_meas_0_series_method = cantor(level);
%
% The time for the series method is less than the time for the integral
% equation method for level<12. For level>=12, the integral equation method
% is faser than the series method.
%
har_meas_0_series_method
har_meas_0_ie_method
%
%