% UNITS
% pressure p in Pa
% volume-vlow v in ???
% length L in m
% diameter d in m
% density rho in ???
% kinematic viscosity nu in ??

b = 6;      % branches
n = 4;      % nodes
c = 3;      % components

B = zeros(b,c);
B_A = [-1,0,-1;0,1,-1;1,1,0] % Äste
B_S = eye(size(B_A)); % Sehnen
B = [B_A,B_S]; % Maschenmatrix
L = [10;10;10]; % length of tubes 
d = [2;2;2]; % diameter of tubes 
rho = [1.5; 1.5; 1.5] % density
%nu = [16;16;16] % (kinematic) viscosity -----> Still needed here 
v_old = [1;1;1;1;1;1]; % Initial vector of volume-flows 

p_mittel = [300000;290000;290000];
T_mittel = [329;344;344];

while true
nu  = CalculateViscosity( p_mittel,T_mittel)
Z = CalculateResistance(b,n,c,B_A,rho,nu,v_old,L,d)

% Creating system of equations
Z_M = B*Z*B';
E_Z = eye(b);
E_K = eye(n-1); % size(E_K)
V_U = zeros(b);
Z_V = zeros(b);
null_one = zeros(b);
null_two = zeros(c);
null_three = zeros(c,b);
middle_A = [E_K,-B_A';null_two,Z_M];
side_A = [null_three;B];
A = [E_Z,-Z,-E_Z;null_one,middle_A,side_A;-V_U,-Z_V,E_Z];

% Creating Solutions-vector
p = zeros(b,1);
p(1) = 400000;
v_rhs_middle = zeros(n-1,1) ;
rhs_middle = [v_rhs_middle;B*p] ;
rhs_end = zeros(b,1);
r_k = [-p;rhs_middle;rhs_end];

% Calculation of the new values 
x = linsolve(A,r_k);

v_new = x(b+1:2*b)

% Check if values are within the tolerance
if  abs(v_old - v_new) <= 0.0001 %perhaps not only compare v, -> whole solution vector
   break;
end
    v_old = v_new
end
