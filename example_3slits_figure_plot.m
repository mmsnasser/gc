clear
% example_2c_figure_plot.m
%
% In this example, plot the filed of the generalized condenser (B,E,delta)
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
n      =  2^11;
t      =  (0:2*pi/n:2*pi-2*pi/n).';
a      =  -0.5; b = 0.5; c = 2;
Lc     =  [-(1+c)/2 ; (a+b)/2 ; (1+c)/2 ];
Lk     =  [ c-1     ;  b-a    ;  c-1    ];
thetk  =  [ 0       ;  0      ;  0      ];
%
[et , etp , cent] = PreImageStrSlit (Lc , Lk , thetk , 0.5 , n , 1e-14 , 100 );
%
%%
figure;
hold on
box on
k=1;  crv = [-c -1];
plot(real(crv),imag(crv),'b','LineWidth',1.2)
k=2;  crv = [a b];
plot(real(crv),imag(crv),'b','LineWidth',1.2)
k=3;  crv = [1 c];
plot(real(crv),imag(crv),'b','LineWidth',1.2)
axis equal
axis([-2.25 2.25 -1 1])
% axis off
set(gca,'XTick',[-2:2],'FontSize',18);
set(gca,'YTick',[-1:1]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_3slits_s
%%
figure;
hold on
box on
k=1;  crv = et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
plot(real(crv),imag(crv),'b','LineWidth',1.2)
k=2;  crv = et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
plot(real(crv),imag(crv),'b','LineWidth',1.2)
k=3;  crv = et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
plot(real(crv),imag(crv),'b','LineWidth',1.2)
axis equal
axis([-2.25 2.25 -1 1])
% axis off
set(gca,'XTick',[-2:2],'FontSize',18);
set(gca,'YTick',[-1:1]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_3slits_e
%%