function [cap , uz] = capgc(et,etp,alphav,deltav,m,mp,ell,alpha,z)
% Compute the capacity of the generalized condensers (B,E,delta)
% 
% Input:
% 1,2) et, etp: parametrization of the boundary and its first derivative
% 3) alphav=[alphav(1),...,alphav(mp)]: alphav(j) is an auxiliary point
% interior to the boundary component \Gamma_j
% 4) deltav=[deltav(1),...,deltav(m)]: deltav(j) is the value of the
% potential function u on \Gamma_j
% 5) m: the number of the closed sets E_k
% 6) mp: mp=m-1 if \Gamma_m is the external boundary component of G,
%    o.w., mp=m 
% 7) ell: the multiplicity of the domain B (B=C for ell=0)
% 8) alpha: for bounded G, alpha is an auxiliary point in G
%           for unbounded G, alpha=inf
% 9) z: a row vector of points in G (if it is required to compute u(z))
% 
% Output:
% cap (the capacity of the generalized condensers (C,E,delta)).
% uz (the values of the potential function u(z) if z is given).
%
% Computing the constants \h_{j,k} for j=1,2,...,m+ell and k=1,2,...,mp
ellp = ell ; ellp(abs(alpha)<inf & mp==m)=ell-1; 
n=length(et)/(m+ell); tht=zeros(size(et)); tht(m*n+1:end)=pi/2;
if mp==m & ellp==ell
    A=exp(-i.*tht);
else
    A=exp(-i.*tht).*(et-alpha);
end
for k=1:mp
    for j=1:m+ell
        jv = 1+(j-1)*n:j*n;
        if (ellp==ell)
            gamk{k}(jv,1)=real(exp(-i.*tht(jv)).*clog(et(jv)-alphav(k)));
        else
            gamk{k}(jv,1)=real(exp(-i.*tht(jv)).*...
                clog((et(jv)-alphav(k))./(et(jv)-alpha)));
        end
    end    
    [mu{k},h{k}]=fbie(et,etp,A,gamk{k},n,5,[],1e-14,100);
    for j=1:m+ell
        hjk(j,k)=mean(h{k}(1+(j-1)*n:j*n));
    end
end
% Computing the constants a_k  for k=1,2,...,m
mat=hjk; mat(1:m,mp+1)=1; mat(m+1:m+ell,mp+1)=0; 
mat(1:m,mp+2:mp+ell+1)=0; mat(m+1:m+ell,mp+2:mp+ell+1)=-eye(ell); 
rhs(1:m,1)=deltav; rhs(m+1:m+ell,1)=0; 
if mp==m
    mat(m+ell+1,1:m)=1; mat(m+ell+1,m+1:m+ell+1)=0; rhs(m+ell+1,1)=0;
end
x=mat\rhs; a=x(1:mp,1); c=x(mp+1);
if mp==m-1
    a(m,1)=-sum(a);
end
% Computing the capacity
cap  = (2*pi)*sum(deltav(:).*a(:));
% compute the values of the potential function u(z) if z is given
if nargin==9
    fet = zeros(size(et)); uz=zeros(size(z));
    for k=1:mp
        fet = fet+a(k).*(gamk{k}+h{k}+i.*mu{k})./A;
        uz=uz-a(k)*log(abs(z-alphav(k)));
    end
    if abs(alpha)<inf
        fz=fcau(et,etp,fet,z);
        uz=uz+c+real((z-alpha).*fz);
    else
        fz=fcau(et,etp,fet,z,n,0);
        uz=uz+c+real(fz);
    end
end
end