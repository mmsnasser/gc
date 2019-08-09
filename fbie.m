function  [mu , h ]  =  fbie(et,etp,A,gam,n,iprec,restart,gmrestol,maxit)
% The function 
%        [mu,h]  =  fbie(et,etp,A,gam,n,iprec,restart,gmrestol,maxit)
% return the unique solution mu of the integral equation 
%               (I-N)mu=-Mgam 
% and the function 
%                h=[(I-N)gam-Mmu]/2,
% where et is the parameterization of the boundary, etp=et', 
% A=exp(-i\thet)(et-alp) for bounded G and by A=exp(-i\thet) for unbounded
% G, gam is a given function, n is the number of nodes in  each boundary
% component, iprec is the FMM precision flag, restart is the maximum number 
% of GMRES method inner iterations, gmrestol is the tolerance of the GMRES   
% method, and maxit is the maximum number of GMRES method outer iterations
%
% Author: Mohamed M S Nasser, v 1.0, 10 December 2017.
% 
% Please cite this function as:
%  M.M.S. Nasser, Fast solution of boundary integral equations with the 
%  generalized Neumann kernel, Electronic Transactions on Numerical 
%  Analysis,  44 (2015) 189--229.
% 
%  PLEASE note that this toolbox contains the files:
%  zfmm2dpart.m
%  fmm2d_r2012a.mexw32
%  fmm2d_r2012a.mexw64
%  pthreadGC2-w32.dll
%  pthreadGC2-w64.dll
%  From the Toolbox:
%  L. G REENGARD AND Z. G IMBUTAS , FMMLIB2D: A MATLAB toolbox for
%  fast multipole method in two dimensions, Version 1.2, 2012.
%  http://www.cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html
%  PLEASE also cite the FMMLIB2D toolbox.
%
%%
a        = [real(et.') ; imag(et.')];
m        =  length(et)/n-1;
b1       = [etp./A].';
[Ub1]    = zfmm2dpart(iprec,(m+1)*n,a,b1,1);
Eone     = (Ub1.pot).';
%%
b(1,1) = 0;
for k=2:n
    b(k,1) = (-1)^(k+1)*(1/n)*cot(pi*(k-1)/n);
end
%%
mu     = gmres(@(x)fB(x),-fC(gam),restart,gmrestol,maxit);
if( nargout == 2 )
    h   = (fC(mu)-fB(gam))./2;
end
%%
function  hx  = fB (x)
    bx2   = [x.*etp./A].';
    [Ubx2]= zfmm2dpart(iprec,(m+1)*n,a,bx2,1);
    Ex    = (Ubx2.pot).';
    hx    =  2.*x-(2/n).*imag(A.*Eone).*x+(2/n).*imag(A.*Ex);
end
%%
function  hx = fC (x)
    bx    = [x.*etp./A].';
    [Ubx] = zfmm2dpart(iprec,(m+1)*n,a,bx,1);
    Ex    = (Ubx.pot).';
    for k=1:m+1
        hLx(1+(k-1)*n:k*n,1) = ifft(fft(b).*fft(x(1+(k-1)*n:k*n,1)));
    end
    hx    = -(2/n).*real(A.*Ex)+(2/n).*real(A.*Eone).*x+hLx;    
end
%%
end