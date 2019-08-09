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
level =  3;
% 
n      =  2^10;
t      = (0:2*pi/n:2*pi-2*pi/n).';
ratio  =  0.5;
%
% 
E     = [0 ; 1];
number_of_intervales = 2^level
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
[x,y] =  meshgrid([-1.0:0.01:2.0],[-1.0:0.01:2.0]);
z     =  x+i.*y;
for k=1:m
    crv    =  et((k-1)*n+1:k*n,1);
    rx=0.5*(max(real(crv))-min(real(crv)));
    ry=0.5*(max(imag(crv))-min(imag(crv)));
    z((real(z-cent(k))./rx).^2+(imag(z-cent(k))./ry).^2<=1)=NaN+i*NaN;
end
[mz,nz]  = size(z); ind = 1;
for k=1:mz
    for j=1:nz
        if abs(z(k,j))>=0
            zv(ind)=z(k,j);
            ind = ind+1;
        end
    end
end
size(zv)
% 
%
[cap,uzv] = capgc (et,etp,alphav,deltav,m,mp,ell,inf,zv);
wzv       = zv+fcau(et,etp,fet,zv,n,0);
%
% 
% 
%%
ind = 1;
for k=1:mz
    for j=1:nz
        if abs(z(k,j))>=0
            uz(k,j)=uzv(ind);
            wz(k,j)=wzv(ind);
            ind = ind+1;
        else
            uz(k,j)=NaN;
            wz(k,j)=NaN;
        end
    end
end
%
%
%%
figure;
hold on
tv = 0:0.05:1;
contour(x,y,uz,tv,'Color','k','LineWidth',1.0,'showtext','on');
for k=1:m+ell
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
axis equal
axis([-1  2 -1  2])
% axis off
box on
ax=gca;
set(gca,'XTick',[-2:2],'FontSize',18);
set(gca,'YTick',[-1:1]);
set(gca,'LooseInset',get(gca,'TightInset'))
%%
% 
figure;
hold on
tv = [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.945855,0.96,0.98];
contour(real(wz),imag(wz),uz,tv,'Color','r','LineWidth',1.0);
for k=1:m
    crv = wet((k-1)*n+1:k*n,1);
    plot(real(crv),imag(crv),'b','LineWidth',2)
end
axis equal
axis([-1  2 -1  2])
% axis off
box on
ax=gca;
set(gca,'XTick',[-2:2],'FontSize',18);
set(gca,'YTick',[-1:1]);
set(gca,'LooseInset',get(gca,'TightInset'))
% print -depsc fig_cant_set_lc
print -dpdf cant_set_3
%