clear
% example_2c_figure_plot.m
%
% In this example, plot the filed of the generalized condenser (B,E,delta)
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
tic
level  =  2;
n      =  2^10;
t      = (0:2*pi/n:2*pi-2*pi/n).';
%
%
ver_out = [1+i ; i ; 0 ; 1];
ver_in  = [1+i ; 1-i ; -1-i ; -1+i]./2;
centa   = [-2-2i ; -2i ; 2-2i ; -2 ; 2; -2+2i ; 2i ; 2+2i]./6;
% 
% outer square
k = 2; cent(k,1)=0.5+0.5i; rad(k,1) = 1;
[et((k-1)*n+1:k*n,1) , etp((k-1)*n+1:k*n,1) ]     =  polygonp(ver_out,n/4);
% 
% 
% first level: one inner square
k = 1; cent(k,1)=0.5+0.5i; rad(k,1) = 1/3;
[et((k-1)*n+1:k*n,1) , etp((k-1)*n+1:k*n,1) ]     =  polygonp(cent(k)+ver_in./3,n/4);
vert{2} = cent(k)+ver_in./3;
% 
% 
for kk=2:level
    if kk==2
        % second level: add 8 more inner squares
        cent(3:3+8-1,1) = cent(2)+centa; 
        for k=3:3+8-1
            [et((k-1)*n+1:k*n,1),etp((k-1)*n+1:k*n,1)] = polygonp(cent(k)+ver_in./(3^2),n/4);
            vert{k} = cent(k)+ver_in./(3^2);
        end
    end
    %
    if kk==3
        % third level: add 64 more inner squares
        for jj=1:8
            cent(11+(jj-1)*8:11+jj*8-1,1) = cent(jj+2)+centa./3; 
            for k=11+(jj-1)*8:11+jj*8-1
                [et((k-1)*n+1:k*n,1),etp((k-1)*n+1:k*n,1)] = polygonp(cent(k)+ver_in./(3^3),n/4);
                vert{k} = cent(k)+ver_in./(3^3);
            end
        end
    end
    if kk==4
        % fourt level: add 512 more inner squares
        for jj=1:64
            cent(75+(jj-1)*8:75+jj*8-1,1) = cent(jj+10)+centa./(3^2); 
            for k=75+(jj-1)*8:75+jj*8-1
                [et((k-1)*n+1:k*n,1),etp((k-1)*n+1:k*n,1)] = polygonp(cent(k)+ver_in./(3^4),n/4);
                vert{k} = cent(k)+ver_in./(3^4);
            end
        end
    end
    if kk==5
        % fifth level: add 4096 more inner squares
        for jj=1:512
            cent(587+(jj-1)*8:587+jj*8-1,1) = cent(jj+74)+centa./(3^3); 
            for k=587+(jj-1)*8:587+jj*8-1
                [et((k-1)*n+1:k*n,1),etp((k-1)*n+1:k*n,1)] = polygonp(cent(k)+ver_in./(3^5),n/4);
                vert{k} = cent(k)+ver_in./(3^5);
            end
        end
    end
end
%%
alpha  = 0.2965+0.3335i;
m      =  2;
mp     =  m-1; 
ell    =  sum(8.^(1:(level-1)));
alphav =  0.5+0.5i;
deltav = [0 ; 1];
%%
[x,y] =  meshgrid(linspace(0.001,0.999,1001));
z     =  x+i.*y;
for k=2:m+ell
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
%%
figure;
hold on
box on
plot(real(z),imag(z),'.r','LineWidth',1.2)
for k=1:m
    crv = et((k-1)*n+1:k*n,1);
    plot(real(crv),imag(crv),'k','LineWidth',1.5)
end
for k=m+1:m+ell
    crv = et((k-1)*n+1:k*n,1);
    plot(real(crv),imag(crv),'-.k','LineWidth',1.5)
end
plot(real(alphav),imag(alphav),'or','LineWidth',1.2)
plot(real(alpha),imag(alpha),'dk','LineWidth',1.2,'MarkerfaceColor','y')
axis equal
%
%%
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
%%
%
figure;
hold on
box on
tv = [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.98];
contour(x,y,uz,tv,'LineWidth',1.0);
colormap jet
colorbar 
caxis([0 1])
for k=1:m
    crv = et((k-1)*n+1:k*n,1);
    plot(real(crv),imag(crv),'k','LineWidth',1.5)
end
for k=m+1:m+ell
    crv = et((k-1)*n+1:k*n,1);
    plot(real(crv),imag(crv),'-.k','LineWidth',1.5)
end
axis equal
axis([-0.01 1.01 -0.01 1.01])
% axis off
set(gca,'XTick',[0:0.25:1.0],'FontSize',14);
set(gca,'YTick',[0:0.25:1.0]);
set(gca,'LooseInset',get(gca,'TightInset'))
print -depsc fig_carpet_lc1_2
print -dpdf  fig_carpet_lc1_2
%%