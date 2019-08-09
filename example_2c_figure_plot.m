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
a=2; r2=0.5;
alphav   =  [0 ; a];
deltav   =  [0 ; 1];
m        =   2; 
mp       =   m; 
ell      =   0;
% 
n        =   2^10;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
% 
k=1; et (1+(k-1)*n:k*n,1) = exp(-i.*t);
k=m; et (1+(k-1)*n:k*n,1) = a+r2.*exp(-i.*t);
figure;
hold on
k=1;
crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
% fill(real(crv),imag(crv),[0.8 0.8 0.8])
plot(real(crv),imag(crv),'b','LineWidth',1.2)
k=2;
crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
% fill(real(crv),imag(crv),[0.8 0.8 0.8])
plot(real(crv),imag(crv),'b','LineWidth',1.2)
axis equal
axis([-1.2 2.7 -1.5 1.5])
axis off
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_BC_E2c_fig
%%