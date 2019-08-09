clear
%%
% In this example, the function 
% cap = capgc (et,etp,alphav,deltav,m,ell,mp,alpha)
% is used to compute the capacity of the generalized condenser (B,E,delta)
% where
% 1)  B=C (the complex plane) and hence ell=0
% 2)  E=E1 U E2 where E1 is the interior of the circle |z|=1 and E2 is the
% interior of the circle |z-a|=r.
% 3)  delta=[0 1];
% 4)  For this case, E is bounded and m'=m=2, G is unbounded
% 
% The exact capacity for this case is known and given by
%  cap = 2/(1+sqrt(1-4*r/((1+a-r)*(a+r-1))))-1;
%
%
n        =   2^10;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
%
%
a        =  2;   
r        =   0.5;
delt2vc  = [0.15:0.15:0.9];
alphav   =  [0 ; a];
m        =   2; 
mp       =   m; 
ell      =   4;
alpha    =  -1+i;
%
%
k=1; % parametrization of \Gamma_1, the internal boundary
et (1+(k-1)*n:k*n,1)   =           exp(-i.*t);
etp(1+(k-1)*n:k*n,1)   =       -i.*exp(-i.*t);
k=2; % parametrization of \Gamma_2, the internal boundary
et (1+(k-1)*n:k*n,1)   =    a+r.*exp(-i.*t);
etp(1+(k-1)*n:k*n,1)   =   -i*r.*exp(-i.*t);
k=3; % parametrization of L_3, the internal boundary
et (1+(k-1)*n:k*n,1)   =   2i+0.9.*exp(-i.*t);
etp(1+(k-1)*n:k*n,1)   =    -0.9i.*exp(-i.*t);
k=4; % parametrization of L_4, the internal boundary
et (1+(k-1)*n:k*n,1)   =  -2i+0.9.*exp(-i.*t);
etp(1+(k-1)*n:k*n,1)   =    -0.9i.*exp(-i.*t);
k=5; % parametrization of L_4, the internal boundary
et (1+(k-1)*n:k*n,1)   =   -2+0.9.*exp(-i.*t);
etp(1+(k-1)*n:k*n,1)   =    -0.9i.*exp(-i.*t);
k=6; % parametrization of L_5, the external boundary
et (1+(k-1)*n:k*n,1)   =        3.*exp(i.*t);
etp(1+(k-1)*n:k*n,1)   =       3i.*exp(i.*t);
%
% 
for itr=1:length(delt2vc)
    deltav   =  [0 ; delt2vc(itr)];
    cap(itr) = capgc (et,etp,alphav,deltav,m,mp,ell,alpha);
    fprintf('The integral int_G|nabla u|^2dm = %12.6f\n',cap(itr))
end
%
format long 
[delt2vc.' cap.']
%
% 