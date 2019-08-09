% example_square_in_square.m
% Nasser, June 9, 2019
clear; 
% To compute the capacity of the square in square domain in
% Section 4.12 of the paper:
% *** 
%
% 
level  =  1;
n      =  2^9;
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
% 
tic
cap = capgc (et,etp,alphav,deltav,m,mp,ell,inf);
fprintf('The integral int_G|nabla u|^2dm = %12.12f\n',cap)
toc
% 
%
%%
% plot the domain
figure
hold on
box on
for k=1:m
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
plot(real(alphav),imag(alphav),'or','LineWidth',1.2)
axis equal
axis([-0.5 1.5 -0.5 1.5])
% axis off
set(gca,'XTick',[-0.5:0.5:1.5],'FontSize',18);
set(gca,'YTick',[-0.5:0.5:1.5]);set(gca,'LooseInset',get(gca,'TightInset'))
%%
% 
% 
capv = [
    4    4.652550988612     1.38
    16   4.562141316258     10.94
    64   4.531268395465     153.24
    256  4.519885907207     2371.08
    ];