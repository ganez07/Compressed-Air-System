
function [Z,z_a] = Copy_of_CalculateResistance(b,n,c,sections,B_A,rho,nu,m_old,L,d,deltaP, P_q, Gesamtknotenindex,p_knoten,m1_real,delta_P_sehne)
% CalculateResistance calculates the resistance of all branches

%, P_q, P_st_old


% INPUT:
%   - branches b
%   - nodes n
%   - loops c
%   - B_A (loop-matrix, Äste)
%   - density rho [kg/m^3]
%   - v_old (inital solutions-vector)
%   - Length of pipe L [m]
%   - diameter of pipe d [m]

% OUTPUT:
%   - resistance matrix Z

B_S = eye(size(B_A));
% Calculating the resistance of all "real" branches (Äste)

max_sections = max(sections); %max. number of sections 

R_e = zeros(n-1,max_sections);  % zeros in the "non-existing" sections
z_a = zeros(n-1,max_sections);  % zeros in the "non-existing" sections
lambda = zeros((n-1),max_sections);

for i=1:(n-1) 
    
for k = 1:sections(i)
    
    R_e(i,k) = (4/(pi*d(i,k)*nu(i,k)))*(abs(m_old(i)/rho(i,k)));  % calcultaion of the Reynols-Number % k = sections, i = nodes
    
    if R_e(i) < 2320                        % calculation of lambda (3 different equations/cases)
        lambda(i,k) = 64/R_e(i,k);    
    elseif 2320 <= R_e(i,k) < 10^5 
        lambda(i,k) = 0.3164*R_e(i,k).^(-0.25);   
    elseif  10^5 <= R_e(i,k) < 10^6
        lambda(i,k) = 0.0032 + (0.221*R_e(i,k).^(-0.237)); 
    else 
        syms lam
        f = 2*log(R_e(i,k)*sqrt(lam))- 0.8 - 1/(sqrt(lam));
        lambda(i,k) = solve(f,lam);
    end 

    
    z_a(i,k) = lambda(i,k) *(8*rho(i,k)/pi.^2)*(L(i,k)/d(i,k).^5)*(m_old(i)/rho(i,k)^2)
end
end

Ges_RES= sum(z_a,2);% calculates the sum of every row -> Gesamtwiderstand der einzelnen Zweige (als Vektor)


% Resistance values for the "imaginary" branches (Sehnenwiderstände Z_S)
Z_A=zeros(n-1); % To calculate the resistance of the "imaginary" branches (Sehnen) a resistance-Matrix of branches (Äste) is needed
for i=1:(n-1)
        Z_A(i,i) = Ges_RES(i);
end
% 
% Zaehler = B_A*Z_A*m_old(1:(n-1))-B_A*P_q(1:(n-1))%+B_A*P_st_old(1:(n-1))%  % fehlen Druckquellen und gesteuerte Druckabfälle: -B_A*P_q(1:(n-1)) +B_A*P_st_old(1:(n-1))
% Nenner = B_S*m_old((n):b)
   Z_s = zeros(c,1); 
% for i=1:c
%        if m_old(i) == 0
% 
%         Z_s(i)= 0
%        else 
%         Z_s(i)= Zaehler(i)/Nenner(i);  
%        end
% end

for i = n:b
     if m_old(i) == 0

        Z_s(i)= 0
     else 
        Z_s(i)=delta_P_sehne(i-n+1)/(m_old(i)) %(m1_real/m_old(1))^(-1)
     end
   
end

  Z_s = Z_s(n:b)

% for i=1:(b-n+1)
% Z_s(i) = abs(delta_P_sehne(i)/(m_old(i+n-1)))
% end

% "Richtiger Massenstrom

%                 p_knoten_theoretisch = p_knoten
% 
%                 for i = 2:n
% 
%                 p_knoten_theoretisch(i) = abs((m1_real/m_old(1))^(-1))*p_knoten(i)
% 
%                 end
% 
% 
%                 Z_s = zeros((b-n+1),1)
%                 deltaP_Sehne = zeros((b-n),1)
% 
%                 for i = n:b
%                 deltaP_Sehne(i-n+1) = Gesamtknotenindex(i,:)*p_knoten_theoretisch %*(-sign(m_old(1)))
%                 
%              
%                 Z_s(i-n+1) = (deltaP_Sehne(i-n+1))/m_old(i)
%             
%                 end




% Adding the resistance of Äste and Sehnen to get a matrix with the resistance of all branches
Z_as = [Ges_RES; Z_s]

Z = zeros(b);

Z_as

for i = 1:b
    Z(i,i) = Z_as(i);
end
end

