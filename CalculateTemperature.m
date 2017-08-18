function [T_node] = CalculateTemperature(T1,n,knotenindex,p_node,sections)
% CalculateTemperature calculates the temperature after every component/pipesection (in Kelvin)
% assumption: equation of state is under reversibel adiabatic condition(same as isentropic case).

% INPUT
%   - pemperature at first node T1 [K]
%   - number_of_nodes n
%   - knotenindex 
%   - pressure at the nodes p_node [Pa]

% OUTPUT
%   - temperatrure at the nodes [K]


T_node = zeros(n-1,max(sections)*2);
T_node (1,1) = T1;      % in K

cp = 1.005;    % specific heat capacity (isobar) for air [kJ/(kg*K)] 
cv = 0.718;    % specific heat capacity (isochor) for air [kJ/(kg*K)] 

n_exponent = cp/cv;  % n_exponent is the necessary coefficient for temperatur's calculation.

%% Calculates the temperature at the node, ignore 'sections'.
    for i = 2:(n-1)
        for k = 1:(n-1)
            if knotenindex((i-1), k) == -1        % to choose the 'real' flow directions, except T1(already exists).
   
                T_node(i,1) = (T_node(k)*(p_node(k,1)/p_node(i,1))^((1-n_exponent)/n_exponent)); % Vorsicht bei der Nummerierung der Ströme und Knoten -> Fehler wenn Strom 3 aus Knoten 4 fließen würde
            
            end
        end  
    end
   
   
%% Calculates the temperature after every component/section, include 'sections' at each branch.
% to calculate the subtemperatur at each 'section' 
    for i=1:(n-1)
        for k= 2:(sections(i)*2)
            T_node(i,k) = T_node(i,k-1)*(p_node(i,k-1)/p_node(i,k))^((1-n_exponent)/n_exponent);
        end
    end
    
end

