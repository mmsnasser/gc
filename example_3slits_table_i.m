clear
% example_2c_error_plot.m
%
% In this example, the function 
% cap = capgc (et,etp,alphav,deltav,m,ell,mp,alpha)
% is used to compute the capacity of the generalized condenser (B,E,delta)
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
t      =  (0:2*pi/n:2*pi-2*pi/n).';
av     =  [-0.9 ; -0.5 ; -0.9 ; 0   ; -0.5 ; -0.7 ; 0.5 ];
bv     =  [ 0   ;  0.5 ;  0.9 ; 0.9 ;  0.5 ;  0.2 ; 0.8 ];
cv     =  [ 2   ;  2   ;  2   ; 2   ;  3   ;  3   ; 3   ];
% 
m        =   2; 
mp       =   m; 
ell      =   1;
deltav   =  [0 ; 1 ];
% 
for k=1:length(av)
    a = av(k); b = bv(k); c = cv(k);
    %
    Lc     =  [-(1+c)/2 ; (1+c)/2 ; (a+b)/2+i ];
    Lk     =  [ c-1     ;  c-1    ;  b-a    ];
    thetk  =  [ 0       ;  0      ;  0      ];    
    %
    [et , etp , cent] = PreImageStrSlit (Lc , Lk , thetk , 0.5 , n , 1e-14 , 100 );
    %
    alphav   =  cent(1:2);
    %
    cap(k,1) = capgc (et,etp,alphav,deltav,m,mp,ell,inf);
end
%
format long 
[av bv cv  cap]
%
%
% 