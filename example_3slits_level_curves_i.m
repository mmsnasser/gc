clear
% example_2c_error_plot.m
%
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
%%
n      =  2^11;
t      =  (0:2*pi/n:2*pi-2*pi/n).';
a      =  -0.5; b = 0.5; c = 2;
Lc     =  [-(1+c)/2 ; (1+c)/2 ; (a+b)/2+i ];
Lk     =  [ c-1     ;  c-1    ;  b-a    ];
thetk  =  [ 0       ;  0      ;  0      ];
%
[et , etp , cent , fet] = PreImageStrSlit (Lc , Lk , thetk , 0.5 , n , 1e-14 , 100 );
%
alphav   =  cent(1:2);
deltav   =  [0 ; 1];
m        =   2; 
mp       =   m; 
ell      =   1;
%
%
[x,y] =  meshgrid([-2.5:0.01:2.5],[-1.1:0.01:2.1]);
z     =  x+i.*y;
for k=1:3
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
cap
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
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
axis equal
axis([-2.5  2.5 -1  2])
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
tv = [0.025,0.05:0.05:0.35,0.405,0.45,0.5,0.55,0.595,0.65:0.05:0.95,0.975];
contour(real(wz),imag(wz),uz,tv,'LineWidth',1.0);
colormap jet
colorbar 
caxis([0 1])
k=1;  crv = [-c -1];
plot(real(crv),imag(crv),'k','LineWidth',1.5)
k=2;  crv = [a b]+i;
plot(real(crv),imag(crv),'-.k','LineWidth',1.5)
k=3;  crv = [1 c];
plot(real(crv),imag(crv),'k','LineWidth',1.5)
axis equal
axis([-2.5  2.5 -1  2])
% axis off
box on
ax=gca;
set(gca,'XTick',[-2:2],'FontSize',14);
set(gca,'YTick',[-1:2]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_3slits_lc_gi
% print -dpdf  fig_3slits_lc_gi
%