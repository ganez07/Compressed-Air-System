function [ p_k,p_knoten ] = CalculateNodepressure(deltaP_ges, p1,p_const, n, knotenindex, sections, z_sections, massenstroeme)
% CalculateNodepressure calculates the pressure at every node.
% The pressure at node 1 (p1) is given (atmosphere)

% INPUT
%   - deltaP (pressure drop)
%   - p1 (pressure at node 1) [Pa]
%   - number of nodes n
%   - knotenindex


%   - deltaP_ges ist nämlich deltaP im Lösungsvektor
%   - deltaP_pipe ist die einzelne Druckabfall in jeweiligen Zweig, "sections considered"
%   - deltaP_matrix ist die Zusammenfassung von deltaP_pipe.



% OUTPUT
%   - pressure at every node [Pa]


p_k = zeros(n-1,max(sections)*2);  % max(sections)*2 to prove, whether there are 'sections' on this branch. If not, zeros in "non-existing" sections

p_k (1,1) = p1;


%% Calculate the pressure at every node, (ignore sections)

p_knoten = zeros(n,1);  
p_knoten (1) = p1;

    for i = 2:n      
        for k = 1:n
            if knotenindex((i-1), k) == -1        % choose the 'real' flow direction.   
   
                p_knoten(i) = abs( p_knoten(k)-deltaP_ges(i-1)); % Vorsicht bei der Nummerierung der Ströme und Knoten -> Fehler wenn Strom 3 aus Knoten 4 fließen würde
            
            end
        end  
    end

%% Chooses the 'real' branch(in which node the flow flows out of the node), (sections considered)
    for i = 2:(n-1)
        for k = 1:(n-1) 
            if knotenindex(i, k) == -1              
                     p_k(i,1) =( p_knoten(k));

            end
        end  
    end
    
%% Calculate the pressure drop of each 'Rohrabschnitt', (sections considered)
deltaP_pipe = zeros(n,max(sections));

for i= 1:n-1
    for k = 1:sections(i)
        deltaP_pipe (i,k) = z_sections(i,k)*massenstroeme(i);
    end
end

%% Matrix with all pressure drops
deltaP_matrix=  zeros(n,max(sections)*2);
for i = 1:(n-1)
 if sections(i)>1
        for q = 2:2:(sections(i)*2-2)
            deltaP_matrix(i,q)=p_const(i,(q/2));    %konstante Quellen: positiv (Kompressor,...); negativ: Verbraucher  
        end  
 end
    
     for k = 1:2:(sections(i)*2-1)
               deltaP_matrix(i,k)=deltaP_pipe(i,(k+1)/2);    % hier Betrag, da durch die Iterationen eventuell ein negativer Knotendruck entstehen könnte was zu imaginären Temperaturen führen würde -> die viskossität wäre dann nicht mehr berechenbar  
     end
       
end

%% Calculates the pressure after every pipesection/component
  for i=1:(n-1)     
       for k = 2:(sections(i)*2)
                p_k(i,k)=abs(p_k(i,k-1)+deltaP_matrix(i,k-1));    % hier Betrag, da durch die Iterationen eventuell ein negativer Knotendruck entstehen könnte was zu imaginären Temperaturen führen würde -> die viskossität wäre dann nicht mehr berechenbar  
       end
  end     
         
end