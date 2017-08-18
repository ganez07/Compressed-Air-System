count = 0
Schritte= [0; 0; 0; 0]


% for the same network, b,n,c stay same at every timestep.
b_values = [6;6;6;6];      % branches (for each timestep)
n_values = [4;4;4;4];      % nodes
c_values = [3;3;3;3];      % loops (Maschen)
sections_values= [2 2 2 2;1 2 2 1;1 1 2 2];       % numer of pipesections in each 'real' branch (Äste:(size(n-1)). row is real branch, column is timestep.


time_steps = 4;        % timesteps
Gesamtknotenindex = [-1 1 0 0; 0 -1 1 0; 0 -1 0 1; -1 0 0 1; 0 0 -1 1; 1 0 -1 0 ]
knotenindex = [-1 1 0 0; 0 -1 1 0; 0 -1 0 1 ];    % row is volume-flow, column is the ith node. it represents the flowdirection. (-1: out of the node, 1: in the node, 0: no relation).
B_A = [-1,0,-1;0,1,-1;1,1,0];                     % loop-matrix for the "real" branches (Äste)

T1_values = [275;300;313;350];                    % temperature at node 1 for every timestep
p1_values = [1*10^5;1*10^5;1*10^5;1*10^5];        % pressure at node 1 for every timestep

m1 = 50
p_knoten = [400; 330; 320; 300]

L1 = [10 10 ;10 0;10 0];  % length of tubes in timestep 1
L2 = [10 10 ;5 5;10 0];  
L3 = [10 10 ;5 5;5 5];  
L4 = [10 10 ;10 0;5 5];  


d1 = [0.2 0.2 ;0.2 0;0.2 0];     % diameter of tubes in timestep 1
d2 = [0.2 0.2 ;0.2 0.2;0.2 0];
d3 = [0.2 0.2 ;0.2 0.2;0.2 0.2];
d4 = [0.2 0.2 ;0.2 0;0.2 0.2];

% Konvergenz überprüfen:



ConvergeCheck = zeros(18,8);

delta_P_sehne_values = [60000, 20000, 30000, 50000; 0, -30000, 0, 30000; 60000, 50000, 30000, 20000]

L_values= [L1,L2,L3,L4]  % to generate matrix of tubelength include every timesteps.

d_values= [d1,d2,d3,d4]  % to generate matrix of tubediameter include every timesteps.

abschnitt= 0; % Variable wird später zur Auswahl des richtigen Teilabschnitts benötigt (Zur Wahl/Abgrenzung der aktuellen Zeitschritte)
%% Initial values                         
x_old = [12;1;1;133;10;10;50;100;2;3;4;-5;1;1;1;1;1;1];  % initial vectors.

%% for-loop to calculate values for each timestep
for t = 1:time_steps        

    x_old = [12;1;1;133;10;10;50;100;2;3;4;-5;1;1;1;1;1;1];
    
    b = b_values(t);                 % branches (to select the  'corresponded vector' for each timestep)
    n = n_values(t);                 % nodes  (for each timestep)
    c = c_values(t);                 % loops  (for each timestep)
    sections= sections_values(:,t)   % sections  (for each timestep)

    
    L=L_values(:,(abschnitt+1):1:(max(sections)+abschnitt)) %wählt den für den aktuellen Zeitschritt gültigen Abschnitt der L_values Matrix aus
    d=d_values(:,(abschnitt+1):1:(max(sections)+abschnitt)) %wählt den für den aktuellen Zeitschritt gültigen Abschnitt der d_values Matrix aus
    
    % Notiz von & für Flo:
    % L_values bzw. L wird durch die Abschnitte immer auf die momentan
    % vorhandenen sections beschränkt. Ästeanzahl bleibt immer gleich.
    % Sections und somit L's sind abhängig von den Prozesssen bzw. der
    % Geometrie des Netzwerkes.
    
    abschnitt = max(sections)+abschnitt % rechnet die Größe des aktuellen Abschitts zur bereits bestehenden Länge -> auch im nächsten Zeitschritt wird wieder der richtige Abschnitt gewählt
    
    delta_P_sehne = delta_P_sehne_values(:,t)
    
    p_values = zeros(b,time_steps);                   % kann sich b mit timesteps ändern???  -> überarbeiten?
    p_values(1,1) = 600000;                           % sources of the branches for every/certain timestep
    p_values(1,2) = 500000;                             
    p_values(1,3) = 500000;
    p_values(1,4) = 500000;
    
    
%     p_values(2,1) = -200000;
%     p_values(3,1) = -200000;
    
    p_values(2,2) = -300000;
    p_values(2,3) = -200000; 
    p_values(3,3) = -200000;
    p_values(3,4) = -300000; 
    
    T_mittel=zeros((n-1),max(sections)) ;        % to create the initial matrix with the right dimension. 
    for i= 1: (n-1)                              % initial assumption: if there are 'sections' in the branch: average temperature T_mittel is 300 [K].
        for k=1:sections(i)                                          % if there are no 'sections' in the branch: T_mittel is 0.
           T_mittel(i,k)=300;
        end
    end

    p_mittel=zeros((n-1),max(sections)) ;        % to create the initial matrix with the right dimension.
    for i= 1: (n-1)                              % initial assumption: if sections exist, p_mittel is 300000 [Pa].if not, 0.
        for k=1:sections(i)
           p_mittel(i,k)=300000;
        end

    end

    rho_x = CalculateDensity(p_mittel,T_mittel,sections);       % density
    T1 = T1_values(t);                                          % temperature at node 1
    %p1 = p1_values(t);                                         % pressure at node 1

    p1=100000
    
    p = p_values(:,t);
    %p(1) = p_values(1,t);                                       % sources of the branches

    %% create important matrices and vectors

    B = zeros(b,c); 
    B_S = eye(size(B_A));           % loop-matrix for the "imaginary" branches (Sehnen)
    B = [B_A,B_S];                  % combined loop-matrix for all braches (Maschenmatrix)
    
for i= 1:5
            ConvergeCheck(:,i) =ConvergeCheck(:,i+1)
            end
            ConvergeCheck(:,6) = x_old
    
    %% main loop to calculate the solutions-vector of the system of equations
        while true
            for i= 1:7
            ConvergeCheck(:,i) =ConvergeCheck(:,i+1);
            end
            ConvergeCheck(:,8) = x_old
            
            m_old = x_old((b+1):(2*b));
            deltaP = x_old(1:(b));
            testbedingung=ones(3*b,1)*0.5
            
           if abs(ConvergeCheck(:,1)+ConvergeCheck(:,3)-ConvergeCheck(:,5)-ConvergeCheck(:,7)) <= 0.0005 && ((abs(ConvergeCheck(:,2)+ConvergeCheck(:,4)-ConvergeCheck(:,6)-ConvergeCheck(:,8))) <= 0.0005) 
% testbedingung && (ConvergeCheck(:,1)+ConvergeCheck(:,3)-ConvergeCheck(:,5)-ConvergeCheck(:,7))> -testbedingung && (ConvergeCheck(:,2)+ConvergeCheck(:,4)-ConvergeCheck(:,6)-ConvergeCheck(:,8)) < testbedingung && (ConvergeCheck(:,2)+ConvergeCheck(:,4)-ConvergeCheck(:,6)-ConvergeCheck(:,8)) > -testbedingung
%            if abs(ConvergeCheck(:,7)/ConvergeCheck(:,5))<0.05&&abs(ConvergeCheck(:,8)/ConvergeCheck(:,6))<0.05
            
           for i = n:b
            m_old (i) = (ConvergeCheck(i,1)+ConvergeCheck(i,3)+ConvergeCheck(i,5)+ConvergeCheck(i,7)+ConvergeCheck(i,2)+ConvergeCheck(i,4)+ConvergeCheck(i,6)+ConvergeCheck(i,8))/8
           
           end
           end
          % 
           
            %% Calculate resistance
            nu = CalculateViscosity(p_mittel,T_mittel,rho_x,sections);
            [Z,z_sections] = Copy_of_CalculateResistance(b,n,c,sections,B_A,rho_x,nu,m_old,L,d,deltaP,p,Gesamtknotenindex,p_knoten,m1,delta_P_sehne);
                            %Copy_of_CalculateResistance(b,n,c,sections,B_A,rho,nu,m_old,L,d,deltaP, P_q, Gesamtknotenindex,p_knoten,m_real)



            %% Create a system of equations
            Z_M = B*Z*B';
            E_Z = eye(b);
            E_K = eye(n-1);    % size(E_K)
            V_U = zeros(b);
            Z_V = zeros(b);
            null_one = zeros(b);
            null_two = zeros(c);
            null_three = zeros(c,b);
            middle_A = [E_K,-B_A';null_two,Z_M];
            side_A = [null_three;B];
            A = [E_Z,-Z,-E_Z;null_one,middle_A,side_A;-V_U,-Z_V,E_Z];  % completed system-matrix

            %% create the solutions-vector
            v_rhs_middle = zeros(n-1,1); 
            rhs_middle = [v_rhs_middle;B*p];
            rhs_end = zeros(b,1);
            r_k = [-p;rhs_middle;rhs_end];    % completed solutinos-vector

            %% Calculate new values to replace initial vector x_old.
            x = linsolve(A,r_k)             % system-solver
            x_new = x; 

            %% Calculate the pressure at all nodes
            pressure_drop= [x_new(1:(n-1))]; 
            mass_flow= [x_new(b+1:(b+n-1))];
            [p_k, p_knoten] = CalculateNodepressure(pressure_drop,p1,p,n,knotenindex,sections,z_sections,mass_flow) %es fehlen noch Variabelen-> deltaP aller Rohrabschnitte -> muss noch programmiert werden
%CalculateNodepressure(deltaP_ges, p1,p_const, number_of_knots, knotenindex, sections, z_sections, mass_flow)
            
            %% Calculating the average pressure in all pipes, only for each section, not for other components:
            % only for each section on each branch
            max_sections = max(sections);           % to choose the largest number of sections.
            p_mittel= zeros(n-1,max_sections);      

            % for i = 1:(n-1)
            % p_mittel(i)= (p_k(i)+p_k(i+1))/2; % wenn wir Druckverluste durch Bauteile in p_Q schreiben, dann: (p_k(1)+p_k(i+1)+p_Q(1)/2)
            % end

            for i = 1:1:(n-1)
                for k = 2:2:(sections(i)*2)                          
                p_mittel(i,(k/2)) =(( p_k(i,k))+p_k(i,(k-1)))/2;
                end
            end


            %% Calculating the temperature at all nodes (sections considered)
            T_node = CalculateTemperature(T1,n,knotenindex,p_k,sections); 

            %% Calculating the average Temperature in all tubes -> for each section!!!
            for i = 1:1:(n-1)
                for k = 2:2:(sections(i)*2)                                 
                   T_mittel(i,(k/2))= (( T_node(i,k))+T_node(i,(k-1)))/2;
                end
            end
            %% Calculating the average density in all tubes
            rho_x = CalculateDensity(p_mittel,T_mittel,sections);

            %% Calculates the volume-flow
            m_real = zeros(size(rho_x));
            for i = 1:n-1
                for k = 1:sections(i)
                    m_real(i,k)= x_old((b+i));
                end
            end
            
            v_real = m_real./rho_x;
           
            count= count+1
            %% Check if values are within the tolerance
            if  abs(x_old - x_new) <= 0.0000000001 
            
                Gesamtlosungsmatrix(:,t) = x_new;   % contains solution-vectors for all timesteps
            %     rho_x % -> Ausgabe für manuelle Berechnung
            %     Z     % -> Ausgabe für manuelle Berechnung
            
            Schritte(t) = count
            count = 0
               break;
            end
                x_old = x_new
                

        end
     
                x_old = x_new
                
 %% hier nochmal die gleiche while-Schleife einfügen mit anderem Calculate resistance, wenn man ein festes m1 vorgegeben hat
end


a = Gesamtlosungsmatrix(:,1)
b = Gesamtlosungsmatrix(:,2)
c = Gesamtlosungsmatrix(:,3)
d = Gesamtlosungsmatrix(:,4)

%% Show results:
    % x_new
    % Druckdifferenzen = x_new(1:b)
    % Massenstroeme = x_new(b+1:2*b)
    % GesteuerteQuellen = x_new(2*b+1:3*b)
    Gesamtlosungsmatrix