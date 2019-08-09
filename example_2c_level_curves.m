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
a        =  2;   
r        =  0.5;
n        =   2^10;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
k=1; % parametrization of \Gamma_1
et (1+(k-1)*n:k*n,1)   =         exp(-i.*t);
etp(1+(k-1)*n:k*n,1)   =     -i.*exp(-i.*t);
k=2; % parametrization of \Gamma_2
et (1+(k-1)*n:k*n,1)   =    a+r.*exp(-i.*t);
etp(1+(k-1)*n:k*n,1)   =   -i*r.*exp(-i.*t);
%
deltav   =  [0 ; 1];
alphav   =  [0 ; a];
m        =   2; 
mp       =   m; 
ell      =   0;
%
%
[x,y] =  meshgrid([-2.5:0.01:4],[-2.5:0.01:2.5]);
z     =  x+i.*y;
z(abs(z)<=1)   = NaN+i*NaN;
z(abs(z-a)<=r) = NaN+i*NaN;
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
%
ind = 1;
for k=1:mz
    for j=1:nz
        if abs(z(k,j))>=0
            uz(k,j)=uzv(ind);
            ind = ind+1;
        else
            uz(k,j)=NaN;
        end
    end
end
%
%
%
figure;
hold on
tv = 0:0.05:1;
contour(x,y,uz,tv,'LineWidth',1.0);
colormap jet
colorbar 
caxis([0 1])
for k=1:m+ell
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',1.5)
end
axis equal
axis([-2.5  4 -2.5  2.5])
% axis off
box on
ax=gca;
set(gca,'XTick',[-2:4],'FontSize',14);
set(gca,'YTick',[-2:2]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_2c_lc
% print -dpdf fig_2c_lc
%
% 