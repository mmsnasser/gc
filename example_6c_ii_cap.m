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
%%
a       =  2;   rv = [0.01,0.05:0.05:0.9,0.925,0.95,0.97,0.98,0.985,0.99];
%% 
n       =   2^10;
t       =  (0:2*pi/n:2*pi-2*pi/n).';
%%
alphav  =  [0 ; a];
deltav  =  [0 ; 1];
m       =   2; 
mp      =   m; 
ell     =   4;
%%
k=1; % parametrization of \Gamma_1, the internal boundary
et (1+(k-1)*n:k*n,1)   =         exp(-i.*t);
etp(1+(k-1)*n:k*n,1)   =     -i.*exp(-i.*t);
k=3; % parametrization of \Gamma_3, the external boundary
et (1+(k-1)*n:k*n,1)   =  1+3i+2.*exp(-i.*t);
etp(1+(k-1)*n:k*n,1)   =     -2i.*exp(-i.*t);
k=4; % parametrization of \Gamma_3, the external boundary
et (1+(k-1)*n:k*n,1)   =  1-3i+2.*exp(-i.*t);
etp(1+(k-1)*n:k*n,1)   =     -2i.*exp(-i.*t);
k=5; % parametrization of \Gamma_3, the external boundary
et (1+(k-1)*n:k*n,1)   =  -2+0.9.*exp(-i.*t);
etp(1+(k-1)*n:k*n,1)   =     -0.9i.*exp(-i.*t);
k=6; % parametrization of \Gamma_3, the external boundary
et (1+(k-1)*n:k*n,1)   =     6+3.*exp(-i.*t);
etp(1+(k-1)*n:k*n,1)   =     -3i.*exp(-i.*t);
% 
for itr=1:length(rv)
    r        =   rv(itr);
    k=2; % parametrization of \Gamma_2, the internal boundary
    et (1+(k-1)*n:k*n,1)   =    a+r.*exp(-i.*t);
    etp(1+(k-1)*n:k*n,1)   =   -i*r.*exp(-i.*t);
    
    % Plot the domain
    figure(1);
    hold on
    box on
    for k=1:m+ell
        crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
        plot(real(crv),imag(crv),'k','LineWidth',1.2)
    end
    axis equal
    drawnow
    
    cap(itr)=capgc (et,etp,alphav,deltav,m,mp,ell,inf);
    fprintf('The integral int_G|nabla u|^2dm = %12.6f\n',cap(itr))    
end
%%
format long g
[rv.' cap.']
%%  