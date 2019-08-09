function [et , etp , cent , fet] = PreImageStrSlit (Lc , Lk , thetk , ratio , n , tolerance , Maxiter )
% PreImageStrSlit.m
% Nasser, June 11, 2019
% Given an unbounded multiply connected domain Omega of connectivity 
% m=length(thetk), border by m segment slits where
% Lc: is a vector of the center points of these slits
% Lc: is a vector of the length of these slits
% thetk: is a vector of the angles between these slits and the positive
% real line.
% 
% This function computes: et , etp , cent , wet; i.e., a preimage 
% unbounded multiply connected domain G border by m ellises:
% et: the parametruzation of the ellipses. 
% etp: the derivative of et.
% cent: the center of the ellipses
% wet: the boundary values of the conformal mapping from G onto Omega,
% 
% ratio, 0<ratio<=1 is the ration of lengthg of minor axix to the major 
% axis of the ellipses. It is better to choose the ratio close to 1. 
% However, when the two slits are close to each other, we need to chose 
% the ratio small but not very small. May be >=0.01. But for small value 
% we need to increase the value of n.
% 
% n: the number of discretization points in each boundary component of G
% tolerance: is the tolerance of the iterative method
% Maxiter: is the maximum number of iterations for the iterative method
%
%
% This code of computing the preimage domain G and the conformal mapping
% is based on the iterative method presented in Section 4 in the paper:
% M. Nasser and C. Green, A fast numerical method for ideal fluid flow 
% in domains with multiple stirrers, Nonlinearity 31 (2018) 815-837.
% 
%%
t         =   (0:2*pi/n:2*pi-2*pi/n).'; 
m         =   length(thetk);
for k=1:m
    thet(1+(k-1)*n:k*n,1)  =  thetk(k);
end
cent     =  Lc;
radx     = (1-0.5*ratio).*Lk;
rady     =  ratio.*radx;
%
err = inf;
itr = 0;
while (err>tolerance)
    itr  =itr+1;  
    for k=1:m
        et(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*thetk(k)).*(+radx(k).*cos(t)-i*rady(k).*sin(t));
        etp(1+(k-1)*n:k*n,1)   =         0.5.*exp(i*thetk(k)).*(-radx(k).*sin(t)-i*rady(k).*cos(t));    
    end
    %
    A       =  exp(i.*(pi/2-thet));
    gam     =  imag(exp(-i.*thet).*et);
    %
    [mun , h ]  =  fbie(et,etp,A,gam,n,5,[],1e-14,200);
    %
    fet            = (gam+h+i.*mun)./A;
    wet              =  et+fet;
    rotwn           =  exp(-i.*thet).*wet;
    for k=1:m
        wnL         =  rotwn((k-1)*n+1:k*n,1);
        centk(k,1)  =  exp(i.*thetk(k)).*((max(real(wnL))+min(real(wnL)))/2+i.*(max(imag(wnL))+min(imag(wnL)))/2);     
        radk(k,1)   =  max(real(wnL))-min(real(wnL)); 
    end
    cent  =  cent -1.0.*(centk-Lc);
    radx  =  radx -(1-0.5*ratio).*(radk-Lk) ;
    rady  =  ratio.*radx;
    err   = (norm(centk-Lc,1)+norm(radk-Lk,1))/m;
    [itr err];
    error (itr,1) = err;
    itrk  (itr,1) = itr;
    %
    if itr>=Maxiter
        'No convergence after Maximunm number of iterations'
        break;
    end
end
%
for k=1:m
    et(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*exp(i*thetk(k)).*(+radx(k).*cos(t)-i*rady(k).*sin(t));
    etp(1+(k-1)*n:k*n,1)   =         0.5.*exp(i*thetk(k)).*(-radx(k).*sin(t)-i*rady(k).*cos(t));    
end
A       =  exp(i.*(pi/2-thet));
gam     =  imag(exp(-i.*thet).*et);
[mun , h ]  =  fbie(et,etp,A,gam,n,5,[],1e-14,200);
fet         = (gam+h+i.*mun)./A;
% wet         =  et+fet;
end