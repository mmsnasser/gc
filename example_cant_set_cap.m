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
for level =  10:10
% 
n      =  2^10;
t      = (0:2*pi/n:2*pi-2*pi/n).';
ratio  =  1;
%
% 
E     = [0 ; 1];
number_of_intervales(level+1,1) = 2^level
for k=1:level
    E  = [E./3 ; E./3+2/3];
end
for k=1:2:length(E)-1
    Lc((k+1)./2,1)=(E(k)+E(k+1))/2;
end
Lk     =  1/(3^level)+zeros(size(Lc));
Lc(end+1,1) = 0.5+i;
Lk(end+1,1) = 1;
thetk  =  zeros(size(Lc));
%
%
%
[et , etp , cent , fet] = PreImageStrSlit (Lc , Lk , thetk , 0.5 , n , 1e-14 , 100 );
wet      =  et+fet;
% 
m        =  length(Lc);
mp       =   m; 
ell      =   0;
alphav   =  cent;
deltav   =  zeros(size(Lc)); deltav(end) = 1;
%
cap(level+1,1)      = capgc (et,etp,alphav,deltav,m,mp,ell,inf)
% 
end
[number_of_intervales cap]
% 
% 
%%
figure;
hold on
box on
for k=1:m
    crv = wet((k-1)*n+1:k*n,1);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
axis equal
axis([-0.5 1.5 -0.5 1.5])
% axis off
set(gca,'XTick',[-2:2],'FontSize',18);
set(gca,'YTick',[-1:1]);
set(gca,'LooseInset',get(gca,'TightInset'))
%%
figure;
hold on
box on
for k=1:m
    crv = et((k-1)*n+1:k*n,1);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
axis equal
axis([-1.00 2.00 -1 2])
% axis off
set(gca,'XTick',[-2:2],'FontSize',18);
set(gca,'YTick',[-1:2]);
set(gca,'LooseInset',get(gca,'TightInset'))
%%