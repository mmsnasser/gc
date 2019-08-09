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
alphav   =  [ 2  ; 2i ; -2  ; -2i ];
deltav   =  [ 1  ; 2  ;  3  ;  4   ; 0];
m        =   5; 
mp       =   m-1; 
ell      =   0;
alpha    =   0;
%
%
% parametrization of \Gamma_1,...,\Gamma_(m-1), the internal boundaries
for k=1:m-1
    et (1+(k-1)*n:k*n,1) = alphav(k)+exp(-i.*t);
    etp(1+(k-1)*n:k*n,1) =    -i.*exp(-i.*t);
end
k=m; % parametrization of \Gamma_m, the external boundary
et (1+(k-1)*n:k*n,1) =   4.*exp(i.*t);
etp(1+(k-1)*n:k*n,1) =  4i.*exp(i.*t);
%
%
[x,y] =  meshgrid([-4:0.01:4],[-4:0.01:4]);
z     =  x+i.*y;
z(abs(z)>=4)   = NaN+i*NaN;
for k=1:m-1
    z(abs(z-alphav(k))<=1) = NaN+i*NaN;
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
[cap,uzv] = capgc (et,etp,alphav,deltav,m,mp,ell,alpha,zv);
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
%%
figure;
hold on
tv = 0.0:0.25:4;
% contour(x,y,uz,tv,'Color','k','showtext','on','LineWidth',1.0);
contour(x,y,uz,tv,'LineWidth',1.0);
colormap jet
colorbar 
caxis([0 4])
for k=1:m+ell
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',1.5)
end
axis equal
axis([-4.1  4.1 -4.1  4.1])
% axis off
box on
set(gca,'XTick',[-4:4],'FontSize',14);
set(gca,'YTick',[-4:4]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_5c_lc
% print -dpdf  fig_5c_lc
%
% 