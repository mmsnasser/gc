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
n        =   2^10;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
%
%
q        =   0.25;
alphav   = [ 0 ];
deltav   = [ 1  ; 0];
m        =   2; 
mp       =   m-1; 
ell      =   0;
alpha    =  (q+1)/2;
%
%
k=1; % parametrization of \Gamma_1, the inner boundary
et (1+(k-1)*n:k*n,1) =    q.*exp(-i.*t);
etp(1+(k-1)*n:k*n,1) = -q*i.*exp(-i.*t);
k=m; % parametrization of \Gamma_m, the outer boundary
et (1+(k-1)*n:k*n,1) =     exp(i.*t);
etp(1+(k-1)*n:k*n,1) =  i.*exp(i.*t);
%
%
%
% [r,thet] =  meshgrid(linspace(0.251,0.999,10),linspace(0,2*pi,40));
[r,thet] =  meshgrid(linspace(0.251,0.999,100),linspace(0,2*pi,200));
x     =  r.*cos(thet); y     =  r.*sin(thet);
z     =  x+i.*y;
[mz,nz]  = size(z); ind = 1;
for k=1:mz
    for j=1:nz
        zv(ind)=z(k,j);
        ind = ind+1;
    end
end
size(zv)
% 
%
[cap,uzv] = capgc (et,etp,alphav,deltav,m,mp,ell,alpha,zv);
cap
%
ind = 1;
for k=1:mz
    for j=1:nz
        uz(k,j)=uzv(ind);
        ind = ind+1;
    end
end
%
uex = log(abs(z))./log(q);
errorz = abs(uz-uex);
norm(errorz,inf)
%
% 
%
%%
figure;
hold on
box on
tv = [0.01,0.05,0.1,0.2,0.5,1,1.2,1.4,1.6,1.8,2].*1e-14;
% contour(x,y,errorz,10,'Color','k','showtext','on','LineWidth',1.0);
contour(x,y,errorz,20,'LineWidth',1.0);
colormap jet
colorbar 
% caxis([1e-16 2e-14])
for k=1:m+ell
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',1.5)
end
axis equal
axis([-1.05  1.05 -1.05  1.05])
% axis off
set(gca,'XTick',[-1:0.5:1],'FontSize',18);
set(gca,'YTick',[-1:0.5:1]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_hm_ann_lc_err1
%%
figure;
hold on
tv = 0:0.1:1;
% contour(x,y,uz,tv,'Color','k','showtext','on','LineWidth',1.0);
contour(x,y,uz,tv,'Color','k','LineWidth',1.0);
for k=1:m+ell
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
axis equal
axis([-1.1  1.1 -1.1  1.1])
% axis off
box on
set(gca,'XTick',[-1:1],'FontSize',18);
set(gca,'YTick',[-1:1]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_hm_ann_lc1
%
%
%%
% figure;
% hold on
% box on
% pcolor(x,y,errorz);
% axis equal
% hold on
% shading interp
% colormap jet
% colorbar 
% set(gca,'XTick',[-1:1],'FontSize',18);
% set(gca,'YTick',[-1:1]);
% set(gca,'LooseInset',get(gca,'TightInset'))
% print -depsc fig_hm_ann_err1
%
%