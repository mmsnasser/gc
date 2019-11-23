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
n      =  2^10;
t      = (0:2*pi/n:2*pi-2*pi/n).';
ratio  =  0.5;
%
% 
level =  4;
number_of_intervales = 2^level
c = 0; r = 1.5; % level of finite Cantor set
for k = 1:level
r = r/3; c = [c+2*r c-2*r];
end
Lc     =  sort(c.');
Lk     =  1/(3^(level-1))+zeros(size(Lc));
thetk  =  zeros(size(Lc));
%
%
%
[et , etp , cent , fet] = PreImageStrSlit (Lc , Lk , thetk , ratio , n , 1e-14 , 100 );
zet  =  et+fet;
%
m        =  length(Lc);
alphav   =  cent;
deltav   =  zeros(size(alphav)); deltav(m/4+1:3*m/4)=1;
mp       =   m; 
ell      =   0;
%%
%
%
tic
[et , etp , cent , fet] = PreImageStrSlit (Lc , Lk , thetk , 1 , n , 1e-14 , 100 );
toc
zet  =  et+fet;
%
%%
figure;
hold on
box on
for k=1:m
    J = (k-1)*n+1:k*n;
    crv = et(J);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
axis equal
%%
figure;
hold on
box on
for k=1:m
    J = (k-1)*n+1:k*n;
    crv = zet(J);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
axis equal
axis equal
axis([-1.55  1.55 -0.25  0.25])
% axis off
box on
ax=gca;
set(gca,'XTick',[-1.5:0.5:1.5],'FontSize',14);
set(gca,'YTick',[-0.25:0.25:0.25]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc  fig_cant_set_har_4
%%
