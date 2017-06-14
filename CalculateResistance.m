function [ Z] = CalculateResistance( b,n,c,B_A,rho,nu,v_old,L,d)
%UNTITLED Calculates the resistance of all branches
%   Detailed explanation goes here

B = zeros(b,c);

B_S = eye(size(B_A));

B = [B_A,B_S]; 

R_e=zeros(n-1,1);
z_i=zeros(n-1,1);

for i=1:(n-1)
R_e(i) = (4/(pi*d(i)*nu(i)))*v_old(i); %klammern;  R_e(i) d(i) nur Äste !!!
    
    if R_e(i) < 2320 
        lambda(i) = 64/R_e(i);
    elseif 2320<R_e(i) < 10.^5 %<= sonst nicht ganzer Bereich abgedeckt
        lambda(i) = 0.3164*R_e(i).^(-0.25);
    else 
        lambda(i) = 0.0032 + (0.221*R_e(i).^(-0.237));
    end
    
       % lambda_old = lambda  % warum?? 
      % for i=1:(n-1)
      z_i(i) = lambda(i) *(8*rho(i)/pi.^2)*(L(i)/d(i).^5)*v_old(i); %L(i) d(i)
      % end
       
end

Z_A=zeros(n-1)
for i=1:(n-1)
        Z_A(i,i) = z_i(i);
        end

%Sehnenwiderstände Z_S
Zaehler = B_A*Z_A*v_old(1:(n-1));
Nenner = B_S*v_old((n):b);

Z_s = zeros(c,1) 
for i=1:c
   Z_s(i)= Zaehler(i)/Nenner(i);
    
end

Z_as =[z_i; Z_s];

for i=1:b
    Z(i,i)=Z_as(i);
end


end

