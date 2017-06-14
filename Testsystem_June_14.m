b = 6;      % branches
n = 4;      % nodes
c = 3;      % components

B = zeros(b,c);

B_A = [-1,0,-1;0,1,-1;1,1,0]

B_S = eye(size(B_A));

B = [B_A,B_S];       % maschenmatrix


lambda_old = 0; %Vektor
lambda = 0; % Vektor, warum unterschied zu lambda_old?
v_old = [5000; 2500; 2500; 1000; 1000; 1000];

L = [10;10;10]
d=[2;2;2]% Durchmesser der Rohre kann sich verändern
rho = [1.5; 1.5; 1.5]% roh ändert sich im System 
nu = [16;16;16] %ebenfalls Vektor

v_old = [500000; 250000; 250000; 1000; 1000; 1000];

while true


Z=CalculateResistance(b,n,c,B_A,rho,nu,v_old,L,d)
Z_M = B*Z*B';

E_Z = eye(b);

E_K = eye(n-1);
% size(E_K)

V_U = zeros(b);
Z_V = zeros(b);

null_one = zeros(b);
null_two = zeros(c);
null_three = zeros(c,b);

middle_A = [E_K,-B_A';null_two,Z_M];
side_A = [null_three;B];

A = [E_Z,-Z,-E_Z;null_one,middle_A,side_A;-V_U,-Z_V,E_Z];

p = zeros(b,1);
p(1) = 400000;
v_rhs_middle = zeros(n-1,1) ;
rhs_middle = [v_rhs_middle;B*p] ;
rhs_end = zeros(b,1);

r_k = [-p;rhs_middle;rhs_end];
        x = linsolve(A,r_k);

    v_new = x(b+1:2*b)

    if  abs(v_old - v_new) <= 0.0001
        break;
    end
    v_old = v_new
    
%     ++counter
end




 

% x = linsolve(A,r_k)

% Z-Value calculation:


% disp (Z)
% disp (v_old)
