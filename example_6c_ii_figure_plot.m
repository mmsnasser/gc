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
a=2; r2=0.5;
%% 
n        =   2^10;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
%%
alphav   =  [0 ; a];
deltav   =  [0 ; 1];
m        =   2; 
mp       =   m; 
ell      =   4;
%%
k=1; % parametrization of \Gamma_1, the internal boundary
et (1+(k-1)*n:k*n,1)   =         exp(-i.*t);
k=2; % parametrization of \Gamma_2, the internal boundary
et (1+(k-1)*n:k*n,1)   =    a+r2.*exp(-i.*t);
k=3; % parametrization of \Gamma_3, the external boundary
et (1+(k-1)*n:k*n,1)   =  1+3i+2.*exp(-i.*t);
k=4; % parametrization of \Gamma_3, the external boundary
et (1+(k-1)*n:k*n,1)   =  1-3i+2.*exp(-i.*t);
k=5; % parametrization of \Gamma_3, the external boundary
et (1+(k-1)*n:k*n,1)   =  -2+0.9.*exp(-i.*t);
k=6; % parametrization of \Gamma_3, the external boundary
et (1+(k-1)*n:k*n,1)   =     6+3.*exp(-i.*t);
% 
% 
figure;
hold on
for k=1:m
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
for k=m+1:m+ell
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'-.r','LineWidth',1.2)
end
axis equal
axis([-4 10 -6 6])
% axis off
box on
ax=gca;
set(gca,'XTick',[-4:2:10],'FontSize',18);
set(gca,'YTick',[-6:2:6]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_6c_ii_fig
%%
