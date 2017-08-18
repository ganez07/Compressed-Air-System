function [rho_x] = CalculateDensity(p_mittel,T_mittel,sections)
% CalculateDensity calculates the density of all pipes
% By calculating compressibility factor of air Z to obtain density rho_x.
% formula of Z is based on virial equations.

% INPUT:
%   - average pressure p_mittel [Pa]
%   - average temperature T_mittel [K]

% OUTPUT:
%   - density rho_x [kg/m^3]


M = 28.949;                       % molar mass [g/mol]
R = 8.3143;                       % molare/universal gas constant [J/mol*K]
rho_x = zeros(size(p_mittel));

% translation of the units (the following equations only work with Celsius & Mpa)
p = p_mittel/(10^6);              % Pa -> Mpa
t = T_mittel-273.15;              % K -> °C


% to calculate density for each section at each branch
for i=1:size(p)
    
    for k= 1:sections(i)          % sections(i) = number of sections at branch "i"

    % calculates Z with P[MPa] and T[°C]
    Z(i,k) = 1.00001 - 5.8057*10^(-3)*p(i,k) + 2.6402*10^(-4)*p(i,k)^2 - 3.3297*10^(-7)*t(i,k) + 1.2420*10^(-4)*p(i,k)*t(i,k) - 2.0158*10^(-6)*p(i,k)^2*t(i,k) + 2.4925*10^(-9)*t(i,k)^2 - 6.2873*10^(-7)*p(i,k)*t(i,k)^2 + 5.4174*10^(-9)*p(i,k)^2*t(i,k)^2;

    % calculates density in g/m^3
    rho_x(i,k) = ((p(i,k)*10^6*M)/(R*(t(i,k)+273.15)*Z(i,k)));

    % calculates density in kg/m^3
    rho_x(i,k) = rho_x(i,k)/1000; 
    
end
end

