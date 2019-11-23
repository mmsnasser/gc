function har_meas_0 = cantor(m) % Green function exterior to a finite approximate Cantor set
%
% We seek the function u(z) that is zero on the real intervals c(j)+r(j)[-1,1]
% j = 1,...,J, and harmonic outside these slits, including at z=infty, except with
% u(z)˜log|z| as z->0. u is expanded as
%
% u(z) = log|z| + C + SUM_{j=1}^J {d(j)*log|wj(z)|
% + SUM_{k=1}^N a(j,k)*real(wj(z)^(-k)),
%
% where wj(z) is a Joukowski map of exterior(slit j) to exterior(unit disk),
% with SUM d(j) = -1; all these coefficients are collected in the vector X.
% The unknowns determined by linear least-squares are C, d(1),...,d(J),
% {a(j,k)}. A has dimensions 1+J*npts by 1+J*(N+1).
% m = 2; 
J=2^m; c = 0; r = 1.5; % level of finite Cantor set
% clc
tic
for i = 1:m
r = r/3; c = [c+2*r c-2*r];
end
N = max(4,6-m); npts = round(1.5*N); % no. expansion terms and sample pts
circ = (1+1e-12)*exp(1i*pi*(.5:npts)'/npts); % roots of unity
z = []; for j = 1:J
z = [z; c(j)+r*(circ+1./circ)/2]; end % sample points on the bndry
A = ones(size(z)); % constant term
for j = 1:J
wj = wz(z,j); A = [A log(abs(wj))]; % logarithmic terms
for k = 1:N, A = [A real(wj.^(-k))]; end % set up least-squares matrix
end
A = [A; zeros(1,1+J*(N+1))];
A(end,2:N+1:end) = 1;
rhs = [-log(abs(z)); -1];
X = A\rhs; % solve least-squares problem
C = X(1); X(1) = []; % extract results
d = X(1:N+1:end); X(1:N+1:end) = []; a = X;

[~,csort]=sort(c.'); %Get the order of c
dd=-d(csort);
har_meas_0 = sum(dd(J/4+1:3*J/4));
time_series_method = toc
% Contour plot
% x = linspace(0,1.8,145); y = linspace(-0.2,0.8,75); [xx,yy] = meshgrid(x,y);
% zz = xx+1i*yy; uu = cantorfun(zz);
% for j = 1:J, slit = c(j)+r*[-1 1]; plot(slit+1e-12i,'-b','linewidth',1.2), hold on, end
% levels = -10.^(-1.2:.1:.2); contour(xx,yy,uu,levels,'k'), axis equal, axis([0 1.8 -.2 .8])
% set(gca,'xtick',0:.4:2,'ytick',-.8:.4:.8,'fontsize',5), op = odeset('events',@event);
% for t = pi*(1:2:23)/48, z0 = .01*exp(1i*t); sol = ode23(@dzdt,[0 100],z0,op);
% z = deval(sol,linspace(0,max(sol.x),300)); plot(z,'k'), plot(conj(z),'k'), end
% plot(0,0,'.r','markersize',16), hold off
% 
function u = cantorfun(z)
u = log(abs(z)) + C;
for j = 1:J
cj = c(j); w = wz(z,j); u = u + d(j)*log(abs(w));
for k = 1:N, wk = w.^(-k); kk = k+(j-1)*N; u = u+a(kk)*real(wk); end
u(abs(w)<=1.01) = NaN;
end
end
function w = wz(z,j)
zc = (z-c(j))/r; sgn = real(zc)>0|(real(zc)==0&imag(zc)>0); sgn = 2*sgn - 1;
w = zc + sgn.*sqrt(zc.^2-1);
end
function g = dzdt(t,z)
g = 1./conj(z);
for j = 1:J
zc = (z-c(j))/r; sgn = real(zc)>0|(real(zc)==0&imag(zc)>0); sgn = 2*sgn - 1;
cw = conj(zc+sgn.*sqrt(zc.^2-1)); dwdzc = conj((1+sgn*zc./sqrt(zc.^2-1))/r);
g = g + d(j)*dwdzc./cw;
for k = 1:N, kk = k+(j-1)*N; g = g - k*a(kk)*dwdzc./cw.^(k+1); end
end
g = g./abs(g);
end
function [val,isterm,dir] = event(t,z)
dir = zeros(J,1); isterm = ones(J,1);
val = zeros(J,1); for j = 1:J; val(j) = abs(wz(z,j))-(1+.01/r); end
end
end