%%% This approach is based on the method used in "Bestimmung der
%%% Kühlgasverteilung in Turbogeneratoren durch Kombination der
%%% Finite-Elemente-Berechnung und der Netzwerkanalyse", Chapter 5 (Zhao)

%%% All Matrices are simular to the ones used in Zhao-Dissertation, it is
%%% highly recommended to have look at the Dissertation while working on
%%% this code.

%%% The network used for this test system is shown in the
%%% Zhao-Dissertation, page 69

%%% Creating specific network-loops
anzahl_zweige = 6;
anzahl_maschen = 3;
anzahl_knoten = 4;
Maschenmatrix_B = [-1,0,-1,1,0,0 ; 0,1,-1,0,1,0 ; 1,1,0,0,0,1;]; %B(i,j): 1:Zweig j gehört zur Masche i und wird positiv durchlaufen. ("1" values that brach j is in loop i (+/- = direction)
    % Maschenmatrix relates to "B" // "B" relates to the "Zweige"
LsgVector = [-200000;0;0;0;0;0;0;0;0;-200000;0;200000;0;0;0;0;0;0]; % Solutions-Vector for the system (vector on the right side)

    Maschenmatrix_BA = zeros(anzahl_maschen, anzahl_knoten-1); % relates to the size of "B(A)" // "B(A)" relates to the "Äste"
    
    for j=1:anzahl_knoten-1 % For-Loop fills "B(A)" with the numbers from "B" (Only "Äste"(lefthalf))
        for i=1:anzahl_maschen
            Maschenmatrix_BA(i,j) = Maschenmatrix_B(i,j);
        end
    end  
    
    Zweigimpedanzmatrix = zeros(anzahl_zweige); % relates to the size of "Z"

    for i=1:anzahl_zweige % For-Loop fills "Z" with resistance-values(from tubes, capacities, etc..)
        Zweigimpedanzmatrix(i,i) = 1;
    end

    Zweigimpedanzmatrix(5,5) = 0;
    
    Maschenimpedanzmatrix = Maschenmatrix_B * Zweigimpedanzmatrix * (Maschenmatrix_B');% relates to "Z(M)"

    Einheitsmatrix_Z = eye(anzahl_zweige); % relates to "E(Z)"
    Einheitsmatrix_K1 = eye(anzahl_knoten-1); % relates to "E(K-1)"
    
    Zv= zeros(anzahl_zweige); %Druckabfall durch Volumenstrom in anderen Zweigen (pressure-drop through the flow in other branches)
    Vu= zeros(anzahl_zweige); %Druckabfall durch Druckabfall in anderen Zweigen (pressure-drop through the pressure-drop in other branches)
   

    
    Nullmatrix_1 = zeros(anzahl_zweige); % Matrix of zeros #1 (dimensions need to fit in the system)
    Nullmatrix_2 = zeros (size(Maschenimpedanzmatrix,1),size(Einheitsmatrix_K1,2)); % Matrix of zeros #2 (dimensions need to fit in the system)
    
    Minnere = [Einheitsmatrix_K1 -Maschenmatrix_BA' ; % 4x4 Hypermatrix inside of the main matrix (consiting out of 4 other matrices)
               Nullmatrix_2 Maschenimpedanzmatrix];
    
    Nullmatrix_3 = zeros ((size(Minnere,1))-size(Maschenmatrix_B,1),size(Einheitsmatrix_Z,2)); % Matrix of zeros #3 (dimensions need to fit in the system)
    
    Minnere_2 = [Nullmatrix_3;Maschenmatrix_B]; % combination of "B" and [0]-Matrix #3
     
    Mges = [Einheitsmatrix_Z -Zweigimpedanzmatrix -Einheitsmatrix_Z; % Main matrix describing the whole system
            Nullmatrix_1 Minnere Minnere_2;
            -Vu -Zv Einheitsmatrix_Z];


             
% % Finding the solution of the system with "fsolve" and genEq.m-File   
options = optimoptions('fsolve','Algorithm','trust-region-dogleg','MaxFunctionEvaluations',10000,'MaxIterations',10000)
N = size(Mges,1);  % Number of the columns of the main matrix, needed for the creation of the start-value-vector      
x0 = ones(N,1); % initial-value-vector
x0 = x0+5000; % Adjustment of the initial-vector for fsolve
eqn = @(x) genEq_mini(x,Mges,anzahl_zweige,LsgVector); % function, needed for fsolve
[x, Fval] = fsolve(eqn,x0,options) % should solve the system




%%% Splitting of the solutions-vector in delta_p, V; p_steuer

delta_p = zeros(anzahl_zweige,1);
for i= 1:length(delta_p)
    
   delta_p(i) = x(i);
end

Volumenstrome = zeros(anzahl_zweige,1);
for i= 1:length(Volumenstrome)
    
   Volumenstrome(i) = x(anzahl_zweige+i);
end

p_steuer = zeros(anzahl_zweige,1);
for i= 1:length(p_steuer)
    
   p_steuer(i) = x(2*anzahl_zweige+i);
end

delta_p
Volumenstrome
p_steuer

%%% TESTPLATZ (folgende 3 Zeilen stellen das gesamte GLS graphisch dar)
% syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18;
% b=[x1; x2; x3; x4; x5; x6; x7^2; x8^2; x9^2; x10^2; x11^2; x12^2; x13; x14; x15; x16; x17; x18];
% Funktionen = (Mges*b)-LsgVector
