% example_square_in_square.m
% Nasser, June 9, 2019
clc;clear; 
% To compute the capacity of the square in square domain in
% Section 4.12 of the paper:
% *** 
%
% 
a    =   0.2;
b    =   0.7;
% choose the value of n
n    =   3*2^13
t    =  (0:2*pi/n:2*pi-2*pi/n).';
% The parametization of the boundary
rec_ver    =  [1+i   ; -1+i   ; -1-i    ;   1-i  ]; % Vertices of the outer square
tri1_ver   =  [0+a*i  ; -(b-a)/sqrt(3)+b*i ;  (b-a)/sqrt(3)+b*i]; % Vertices of the first triangle
tri2_ver   =  [0-a*i  ;  (b-a)/sqrt(3)-b*i ; -(b-a)/sqrt(3)-b*i]; % Vertices of the second triangle
[et(1:n,1)    ,etp(1:n,1)    ]     =  polygonp(tri1_ver,n/3);
[et(n+1:2*n,1),etp(n+1:2*n,1)]     =  polygonp(tri2_ver,n/3);
[et(2*n+1:3*n,1),etp(2*n+1:3*n,1)] =  polygonp(rec_ver,n/4);
%
%
alphav   =  [ mean(tri1_ver) ; mean(tri2_ver) ];
deltav   =  [ 1  ; 1  ; 0  ];
m        =   3; 
mp       =   m-1; 
ell      =   0;
alpha    =   0;
% 
% 
[x,y] =  meshgrid([-0.99:0.01:0.99],[-0.99:0.01:0.99]);
z     =  x+i.*y;
[in1 on1] = inpolygon(x,y,real(tri1_ver),imag(tri1_ver));
z(in1)=NaN+i*NaN; 
z(on1)=NaN+i*NaN;
[in2 on2] = inpolygon(x,y,real(tri2_ver),imag(tri2_ver));
z(in2)=NaN+i*NaN; 
z(on2)=NaN+i*NaN;
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
%%
figure;
hold on
tv = 1-[0.08,0.14575,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.96];
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
print -depsc fig_2tri_squ_lc
% print -dpdf fig_2tri_squ_lc
%%
% 
% 