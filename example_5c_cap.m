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
n        =   2^10
t        =  (0:2*pi/n:2*pi-2*pi/n).';
%
%
alphav   =  [ 2  ; 2i ; -2  ; -2i ];
deltav   =  [ 1  ; 2  ;  3  ;  4   ; 0];
m        =   5; 
mp       =   m-1; 
ell      =   0;
alpha    =   0;
%
%
% parametrization of \Gamma_1,...,\Gamma_(m-1), the internal boundaries
for k=1:m-1
    et (1+(k-1)*n:k*n,1) = alphav(k)+exp(-i.*t);
    etp(1+(k-1)*n:k*n,1) =    -i.*exp(-i.*t);
end
k=m; % parametrization of \Gamma_m, the external boundary
et (1+(k-1)*n:k*n,1) =   4.*exp(i.*t);
etp(1+(k-1)*n:k*n,1) =  4i.*exp(i.*t);
%
%
cap = capgc (et,etp,alphav,deltav,m,mp,ell,alpha);
format long
cap
% 
%
figure;
hold on
for k=1:m
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
% plot(real(alpha),imag(alpha),'ok','MarkerfaceColor','k')
axis equal
axis([-4.1  4.1 -4.1  4.1])
% axis off
box on
ax=gca;
set(gca,'XTick',[-4:4],'FontSize',18);
set(gca,'YTick',[-4:4]);
set(gca,'LooseInset',get(gca,'TightInset'))
%
% 