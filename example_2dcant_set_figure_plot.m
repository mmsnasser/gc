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
n      =  2^11;
t      = (0:2*pi/n:2*pi-2*pi/n).';
ratio  =  0.5;
%
% 
E     = [0 ; 1];
level =  3;
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
k = 4^level+1;
et((k-1)*n+1:k*n,1)  = 0.5+0.5i+exp(i.*t);
etp((k-1)*n+1:k*n,1) =       i.*exp(i.*t);
%
figure;
hold on
box on
for k=1:4^level+1
    crv = et((k-1)*n+1:k*n,1);
    plot(real(crv),imag(crv),'b','LineWidth',1.2)
end
axis equal
axis([-0.55 1.55 -0.55 1.55])
% axis off
set(gca,'XTick',[-0.5:0.5:1.5],'FontSize',18);
set(gca,'YTick',[-0.5:0.5:1.5]);
set(gca,'LooseInset',get(gca,'TightInset'))
%
%