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
cap = capgc (et,etp,alphav,deltav,m,mp,ell,alpha);
fprintf('The integral int_G|nabla u|^2dm = %12.8f\n',cap)
% 
%
%
% plot the domain
figure
hold on
for k=1:3
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
plot(real(alphav),imag(alphav),'or','LineWidth',1.2)
axis equal
axis([-1.1 1.1 -1.1 1.1])
% axis off
box on
set(gca,'XTick',[-1:1],'FontSize',18);
set(gca,'YTick',[-1:1]);
set(gca,'LooseInset',get(gca,'TightInset'))
%