% example_square_in_square.m
% Nasser, June 9, 2019
clc;clear; 
% To compute the capacity of the square in square domain in
% Section 4.12 of the paper:
% *** 
%
% 
av     =  [0.1 ; 0.2 ; 0.2 ; 0.3 ; 0.3 ];
bv     =  [0.3 ; 0.4 ; 0.7 ; 0.8 ; 0.9 ];
deltav =  [ 0  ; 0  ; 1  ];
% 
for j=1:length(av)
    a    =   av(j);
    b    =   bv(j);
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
    alphav   =  [ mean(tri1_ver) ; mean(tri2_ver) ];
    m        =   3; 
    mp       =   m-1; 
    ell      =   0;
    alpha    =   0;
    % 
    cap(j,1) = capgc (et,etp,alphav,deltav,m,mp,ell,alpha);
    % 
end
%
format long g
[av bv cap]
%