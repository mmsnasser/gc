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
%%
a        =   2;    r  =  0.5;
n        =   2^10;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
alphav   =  [0  ;  a];
centv    =  [0  ; a  ; 1+3i ; 1-3i ; -2   ;  6  ];
radv     =  [1  ; r  ; 2    ; 2    ;  0.9 ;  3];
deltav   =  [0 ; 1];
m        =   2; 
mp       =   m; 
ell      =   4;
alpha    =  -1+i;
% 
% 
% parametrization of \Gamma_1,...,\Gamma_m+ell, the internal boundaries
for k=1:m+ell
    et (1+(k-1)*n:k*n,1) = centv(k)+radv(k).*exp(-i.*t);
    etp(1+(k-1)*n:k*n,1) =       -i*radv(k).*exp(-i.*t);
end
% 
% 
[x,y] =  meshgrid([-5:0.01:11],[-7:0.01:7.4]);
z     =  x+i.*y;
for k=1:m+ell
    z(abs(z-centv(k))<=radv(k)) = NaN+i*NaN;
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
%%
figure;
hold on
tv = [0.005,0.05,0.1,0.15,0.2,0.23,0.27,0.3,0.33,0.36,0.375,0.383,0.4,...
      0.43,0.47,0.55,0.65,0.75,0.85,0.95,0.97];
contour(x,y,uz,tv,'LineWidth',1.0);
colormap jet
colorbar 
caxis([0 1])
for k=1:m
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',1.5)
end
for k=m+1:m+ell
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'-.k','LineWidth',1.5)
end
axis equal
axis([-5  11 -7  7.2])
% axis off
box on
ax=gca;
set(gca,'XTick',[-5:2:11],'FontSize',14);
set(gca,'YTick',[-7:2:7]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_6c_ii_lc
% print -dpdf  fig_6c_ii_lc
%
% 
%%