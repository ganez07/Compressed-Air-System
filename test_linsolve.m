b = 6;      % branches
n = 4;      % nodes
c = 3;      % components

B = zeros(b,c);

B_A = [-1,0,-1;0,1,-1;1,1,0]

B_S = eye(size(B_A))

B = [B_A,B_S]       % maschenmatrix

Z = zeros (b)       % Zweigimpedanzmatrix

for i=1:b
    Z(i,i) = 1
end

Z_M = B*Z*B'

E_Z = eye(b)

E_K = eye(n-1)
size(E_K)

V_U = zeros(b);
Z_V = zeros(b);

null_one = zeros(b);
null_two = zeros(c);
null_three = zeros(c,b);

middle_A = [E_K,-B_A';null_two,Z_M]
side_A = [null_three;B]

A = [E_Z,-Z,-E_Z;null_one,middle_A,side_A;-V_U,-Z_V,E_Z];

p = zeros(b,1)
p(1) = 20000
rhs_middle = [0;0;0;B*p]
rhs_end = zeros(b,1)

r_k = [p;rhs_middle;rhs_end]

x = linsolve(A,r_k)