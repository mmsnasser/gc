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
a        =   2;    r  =  0.5;
n        =   2^10;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
alphav   =  [0 ; a];
centv    =  [0 ; a ; 2i  ; -2i ; -2  ];
radv     =  [1 ; r ; 0.9 ;  0.9;  0.9];
deltav   =  [0 ; 1];
m        =   2; 
mp       =   m; 
ell      =   4;
alpha    =  -1+i;
%
%
% parametrization of \Gamma_1,...,\Gamma_(m-1), the internal boundaries
for k=1:m+ell-1
    et (1+(k-1)*n:k*n,1) = centv(k)+radv(k).*exp(-i.*t);
    etp(1+(k-1)*n:k*n,1) =       -i*radv(k).*exp(-i.*t);
end
k=m+ell; % parametrization of \Gamma_m, the external boundary
et (1+(k-1)*n:k*n,1) =   3.*exp(i.*t);
etp(1+(k-1)*n:k*n,1) =  3i.*exp(i.*t);
%
% 
% 
% Plot the domain
figure(1);
hold on
box on
for k=1:m+ell
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',1.2)
end
plot(real(alpha),imag(alpha),'or')
axis equal
drawnow
%
% 
cap  =capgc (et,etp,alphav,deltav,m,mp,ell,alpha);
fprintf('The integral int_G|nabla u|^2dm = %12.6f\n',cap)    
%%
