% example_square_in_square.m
% Nasser, June 9, 2019
clear; 
% To compute the capacity of the square in square domain in
% Section 4.12 of the paper:
% *** 
%
% 
level  =  0; 
r      =  128;
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
    [et((k-1)*n+1:k*n,1),etp((k-1)*n+1:k*n,1)] = polygonp(vert{k},n/4);
end
%
k = 4^level+1;
et((k-1)*n+1:k*n,1)  = 0.5+0.5i+r.*exp(i.*t);
etp((k-1)*n+1:k*n,1) =       i.*r.*exp(i.*t);
%
%
%
for k=1:4^level
    alphav(k,1)  =  mean(vert{k});
end
deltav   =  zeros(size(alphav)); deltav(end+1) = 1;
m        =   4^level+1; 
mp       =   m-1; 
ell      =   0;
alpha    =   0.5-0.25i;
% 
% 
cap = capgc (et,etp,alphav,deltav,m,mp,ell,alpha);
fprintf('The integral int_G|nabla u|^2dm = %12.12f\n',cap)
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
plot(real(alpha),imag(alpha),'dr','LineWidth',1.2)
axis equal
% axis off
set(gca,'XTick',[-0.5:0.5:1.5],'FontSize',18);
set(gca,'YTick',[-0.5:0.5:1.5]);set(gca,'LooseInset',get(gca,'TightInset'))
%%
% 
% 
capv = [
    1   11.953076429137    
    2   5.148104923657
    4   3.283370847376
    8   2.410318948828
    16  1.904034193007
    32  1.573518434273
    64  1.340776875382
    128 1.168014094121
    
    4    11.598544818484
    16   11.460681553205
    64   11.408998984224
    256  11.389646462346];