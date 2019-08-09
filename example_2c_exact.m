clear
% example_2c_error_plot.m
%
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
a        =  2;   rv = [0.01,0.05:0.05:0.9,0.925,0.95,0.97,0.98,0.985,0.99];
%
% 
n        =   2^10;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
%
%
alphav   =  [0 ; a];
deltav   =  [0 ; 1];
m        =   2; 
mp       =   m; 
ell      =   0;
%
%
for itr=1:length(rv)
    r        =   rv(itr);
    k=1;
    et (1+(k-1)*n:k*n,1)   =         exp(-i.*t);
    etp(1+(k-1)*n:k*n,1)   =     -i.*exp(-i.*t);
    k=m; % parametrization of \Gamma_m, the external boundary
    et (1+(k-1)*n:k*n,1)   =    a+r.*exp(-i.*t);
    etp(1+(k-1)*n:k*n,1)   =   -i*r.*exp(-i.*t);
    
    % Plot the domain
    figure(1);
    hold on
    box on
    k=1;
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',1.2)
    k=2;
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
    axis equal
    drawnow
   
    exact_q     =  2/(1+sqrt(1-4*r/((1+a-r)*(a+r-1))))-1;
    exact_cap(itr)   =  2*pi/log(1/exact_q)
end
%
%
format long g
[rv.' exact_cap.']
%%
