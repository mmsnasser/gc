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
n      =  2^9;
t      = (0:2*pi/n:2*pi-2*pi/n).';
ratio  =  0.5;
%
% 
level =  6;
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
%
%
%
[x,y] =  meshgrid(linspace(-2,2,1001),linspace(-1,1,1001));
z     =  x+i.*y;
for k=1:m
    crv    =  et((k-1)*n+1:k*n,1);
    rx=0.5*(max(real(crv))-min(real(crv)));
    ry=0.5*(max(imag(crv))-min(imag(crv)));
    z((real(z-cent(k))./rx).^2+(imag(z-cent(k))./ry).^2<=1-1e-3)=NaN+i*NaN;
end
% 
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
    plot(real(crv),imag(crv),'b','LineWidth',1.5)
end
axis equal
axis([-2  2 -1  1])
% axis off
box on
ax=gca;
% set(gca,'XTick',[0:0.4:2],'FontSize',18);
% set(gca,'YTick',[0:0.4:0.8]);
set(gca,'LooseInset',get(gca,'TightInset'))
%%
% 
figure;
hold on
set(gca,'LooseInset',get(gca,'TightInset'))
tv = [0.05:0.05:0.7,0.725,0.75:0.05:0.95];
contour(real(wz),imag(wz),uz,tv,'LineWidth',1.0);
colormap jet
colorbar 
caxis([0 1])
for k=1:m+ell
    crv    =  zet((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',1.5)
end
axis equal
axis([-2  2 -1  1])
% axis off
box on
set(gca,'XTick',[-2:1:2],'FontSize',14);
set(gca,'YTick',[-1:0.5:1]);
set(gcf,'Renderer','opengl')
print -depsc -r1000  fig_cant_set_har_lc6.eps
%%