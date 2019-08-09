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
a        =   2;    r  =  0.5;
n        =   2^10;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
alphav   =  [0 ; a];
centv    =  [0 ; a ; 2i  ; -2i ; -2  ];
radv     =  [1 ; r ; 0.9 ;  0.9;  0.9];
deltav   =  [0 ; 1];
m        =   2; 
mp       =   m; 
ell      =   4;
alpha    =  -1+i;
%
%
% parametrization of \Gamma_1,...,\Gamma_(m-1), the internal boundaries
for k=1:m+ell-1
    et (1+(k-1)*n:k*n,1) = centv(k)+radv(k).*exp(-i.*t);
    etp(1+(k-1)*n:k*n,1) =       -i*radv(k).*exp(-i.*t);
end
k=m+ell; % parametrization of \Gamma_m, the external boundary
et (1+(k-1)*n:k*n,1) =   3.*exp(i.*t);
etp(1+(k-1)*n:k*n,1) =  3i.*exp(i.*t);
%
% 
% 
[x,y] =  meshgrid([-3:0.01:3],[-3:0.01:3]);
z     =  x+i.*y;
z(abs(z)>=3)   = NaN+i*NaN;
for k=1:m+ell-1
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
tv = [0.0005,0.005,0.012,0.02,0.028,0.03335,0.04,0.05,0.1,0.2,0.3,0.4,0.475,...
      0.51,0.59,0.7,0.85,0.95,0.97];
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
axis([-3.1  3.1 -3.1  3.1])
% axis off
box on
ax=gca;
set(gca,'XTick',[-3:3],'FontSize',14);
set(gca,'YTick',[-3:3]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_6c_i_lc
%
% 
%%
