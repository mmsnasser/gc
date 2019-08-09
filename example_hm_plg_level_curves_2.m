% example_square_in_square.m
% Nasser, June 9, 2019
clc;clear; 
% To compute the capacity of the square in square domain in
% Section 4.12 of the paper:
% *** 
%
% 
n        =   5*2^8;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
%
%
et(1:n,1)  =  0.5+0.25.*exp(-i.*t);
etp(1:n,1) = -0.25i.*exp(-i.*t);
et(n+1:2*n,1)  = -0.5+0.25.*exp(-i.*t);
etp(n+1:2*n,1) = -0.25i.*exp(-i.*t);
%
out_ver    =  [1      ;  i     ; -1      ;  -1-i    ;   1-i  ];      % Vertices of the outer ploygon
rec_ver    =  [0.5-0.5i ; 0.5-0.8i ; -0.5-0.8i ; -0.5-0.5i]; % Vertices of the rectangle
[et(2*n+1:3*n,1),etp(2*n+1:3*n,1)] =  polygonp(rec_ver,n/length(rec_ver));
[et(3*n+1:4*n,1),etp(3*n+1:4*n,1)] =  polygonp(out_ver,n/length(out_ver));
%
alphav   = [ 0.5 ; -0.5 ; -0.6i ];
m        =   4; 
mp       =   m-1; 
ell      =   0;
alpha    =   0;
%
deltav   =  [ 0  ; 1  ; 0 ; 0 ];
% 
% 
[x,y] =  meshgrid(linspace(-0.999,0.999,1000),linspace(-0.999,0.999,1000));
z     =  x+i.*y;
[in1 on1] = inpolygon(x,y,real(out_ver),imag(out_ver));
z(~in1)=NaN+i*NaN; 
z(on1)  =NaN+i*NaN;
[in2 on2] = inpolygon(x,y,real(rec_ver),imag(rec_ver));
z(in2)=NaN+i*NaN; 
z(on2)=NaN+i*NaN;
% 
z(abs(z-0.5)<=0.2501) = NaN+i*NaN;
z(abs(z+0.5)<=0.2501) = NaN+i*NaN;
%
% figure;
% hold on
% for k=1:m
%     crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
%     plot(real(crv),imag(crv),'b','LineWidth',1.2)
% end
% plot(real(z),imag(z),'or')
% axis equal
% axis([-1.1  1.1 -1.1  1.1])
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
[~,uzv] = capgc (et,etp,alphav,deltav,m,mp,ell,alpha,zv);
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
%%
figure;
hold on
tv = [0.0000003,0.00001,0.0001,0.0005,0.003,0.01025,0.02,0.05,0.13,0.25,0.4,0.5,0.6,0.7,0.8,0.9,0.98];
contour(x,y,uz,tv,'LineWidth',1.0);
colormap jet
colorbar 
caxis([0 1])
for k=1:m+ell
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',1.5)
end
axis equal
axis([-1.05 1.05 -1.05 1.05])
% axis off
box on
set(gca,'XTick',[-1:0.5:1],'FontSize',14);
set(gca,'YTick',[-1:0.5:1]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_hm_plg_lc_2
% print -dpdf fig_hm_plg_lc_2
%%
% 
% 