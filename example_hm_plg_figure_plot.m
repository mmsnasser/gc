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
n        =   5*2^8;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
%
%
et(1:n,1)  =  0.5+0.25.*exp(-i.*t);
etp(1:n,1) = -0.25i.*exp(-i.*t);
et(n+1:2*n,1)  = -0.5+0.25.*exp(-i.*t);
etp(n+1:2*n,1) = -0.25i.*exp(-i.*t);
%
out_ver    =  [1      ;  i     ; -1      ;  -1-i    ;   1-i  ];      % Vertices of the outer ploygon
rec_ver    =  [0.5-0.5i ; 0.5-0.8i ; -0.5-0.8i ; -0.5-0.5i]; % Vertices of the rectangle
[et(2*n+1:3*n,1),etp(2*n+1:3*n,1)] =  polygonp(rec_ver,n/length(rec_ver));
[et(3*n+1:4*n,1),etp(3*n+1:4*n,1)] =  polygonp(out_ver,n/length(out_ver));
%
%
alphav   = [ 0.5 ; -0.5 ; -0.6i ];
m        =   4; 
mp       =   m-1; 
ell      =   0;
alpha    =   0;
%
%
%
%
figure;
hold on
for k=1:m
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
plot(real(alpha),imag(alpha),'ok','MarkerfaceColor','k')
plot(real(alphav),imag(alphav),'dr','MarkerfaceColor','r')
axis equal
axis([-1.1  1.1 -1.1  1.1])
% axis off
box on
ax=gca;
set(gca,'XTick',[-4:4],'FontSize',18);
set(gca,'YTick',[-4:4]);
set(gca,'LooseInset',get(gca,'TightInset'))
%
% 