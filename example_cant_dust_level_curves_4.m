% example_square_in_square.m
% Nasser, June 9, 2019
clear; 
% To compute the capacity of the square in square domain in
% Section 4.12 of the paper:
% *** 
%
% 
level  =  1;
n      =  2^10;
t      = (0:2*pi/n:2*pi-2*pi/n).';
ratio  =  0.5;
%
% 
E     = [0 ; 1];
number_of_squares= 4^level
for k=1:level
    E  = [E./3 ; E./3+2/3];
end
for k=1:length(E)
    for j=1:length(E)
        F(k,j)=E(j)+i*E(k);
    end
end
ind = 1;
for k=1:2:length(E)
    for j=1:2:length(E)
        vert{ind} = [F(k,j) ; F(k+1,j) ; F(k+1,j+1) ; F(k,j+1)];
        ind = ind+1;
    end
end
%
for k=1:4^level
    [et((k-1)*n+1:k*n,1) , etp((k-1)*n+1:k*n,1) ]     =  polygonp(vert{k},n/4);
end
%
%
%
%
m        =   4^level; 
mp       =   m; 
ell      =   0;
for k=1:m
    alphav(k,1)  =  mean(vert{k});
end
deltav   =   zeros(size(alphav)); 
deltav(m/2+1:m) = 1;
% 
% 
[x,y] =  meshgrid(linspace(-0.5,1.5,1001));
z     =  x+i.*y;
for k=1:4^level
    [in on] = inpolygon(x,y,real(vert{k}),imag(vert{k}));
    z(in)=NaN+i*NaN;
    z(on)=NaN+i*NaN;
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
% plot the domain
figure
hold on
box on
plot(real(z),imag(z),'.r','LineWidth',1.2)
for k=1:m
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
plot(real(alphav),imag(alphav),'or','LineWidth',1.2)
axis equal
axis([-0.55 1.55 -0.55 1.55])
% axis off
set(gca,'XTick',[-0.5:0.5:1.5],'FontSize',18);
set(gca,'YTick',[-0.5:0.5:1.5]);
set(gca,'LooseInset',get(gca,'TightInset'))
% 
% 
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
%%
% plot the domain
figure;
hold on
box on
tv = [0.02,0.032841,0.06,0.1:0.05:0.9,0.94,0.96716,0.98];
contour(x,y,uz,tv,'LineWidth',1.0);
colormap jet
colorbar 
caxis([0 1])
for k=1:m+ell
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',1.5)
end
axis equal
axis([-0.51 1.51 -0.51 1.51])
% axis off
set(gca,'XTick',[-0.5:0.5:1.5],'FontSize',14);
set(gca,'YTick',[-0.5:0.5:1.5]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_cant_dust_lc4
% print -dpdf fig_cant_dust_lc4
%%