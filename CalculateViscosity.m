function [nu] = CalculateViscosity(p_mittel,T_mittel,rho_x,sections)
% CalculateViscosity calculates the kinematic viscosity of each pipe.
% Two different methods have been used under different conditions.

% 1.(small pressure): when average pressure is smaller than 50-100 bar.
% Using the founded formula to calculate kinematic viscosity.
% and A is the important coefficient of the formula.

% 2.(large pressure): when 50 bar < p < 1000 bar. 
% Inperpolation of the values based on the chosen VDI-Wärmeatlas tables.

% Choosing which method depends on different input pressures.
% Result: the input pressure must between 1-1000 bar.

% INPUT: 
%   - average pressure p_mittel [Pa]
%   - average temperature T_mittel [K]
%   - density rho_x [kg/m^3]

% OUTPUT: 
%   - kinematic viscosity kin_vis [m^2/s]*10^-6

% Tranlation of units.
p = p_mittel/(10^5);        % pressure unit: from pascal to bar 
T = T_mittel-273.15;        % temperature unit: from Kelvin to Celsius 

% Initial vector for method 1: formula-calculations
nu = zeros(size(p));

p_c = 39.5;      % critical pressure [bar] (always constant)
T_c = 133;       % critical temperature [K] (always constant)
M = 28.949;      % molar mass [g/mol] (always constant)

vis_1bar = zeros(size(p));
p_r = zeros(size(p));
A = zeros(size(p));
dyn_vis = zeros(size(p));
kin_vis = zeros(size(p));

for i = 1:size(sections)      % size(sections) = size(branches)
    for k = 1:sections(i)     % sections(i) = number of sections at branch "i".
    if p(i,k)>= 1 && p(i,k)<50
    
        % to calculate the dynamic viscosity at p=1 bar in 10^-7 Pas at certain temperature
        % interpolation of the temperature at p=1 bar based on VDI-Wärmeatlas D2.2 Tabelle 13.
        eta=[8.664 10.26 11.78 13.23 14.61 15.94 17.22 18.45 19.64 20.78 21.90 22.98 24.03 25.05 26.05 27.97 29.81 33.28 36.53 39.60 42.52 45.32 48.02 50.63];
        TKord=[-150	-125	-100	-75	-50	-25	0	25	50	75	100	125	150	175	200	250	300	400	500	600	700	800	900	1000];
        F = griddedInterpolant(TKord,eta,'spline');
        vis_1bar(i,k) = (F(T(i,k))*10^(-6))/10^(-7); 
        
        % Calculation of p_r 
        p_r(i,k) = rho_x(i,k)/313;  % the reduced parameter: p_r = rho_x/rho_crit [Kg/m^3]

        % Calcultation of A 
        A(i,k) = 1.023 + 0.23364*p_r(i,k) + 0.58533*p_r(i,k)^2 - 0.40758*p_r(i,k)^3 + 0.093324*p_r(i,k)^4;

        % Main Formula
        dyn_vis(i,k) = vis_1bar(i,k) + (A(i,k)^4 - 1)*((p_c^(2/3)*M^(0.5))/T_c^(1/6));
        dyn_vis(i,k) = dyn_vis(i,k)/10;
        kin_vis(i,k) = dyn_vis(i,k)/rho_x(i,k);
        nu(i,k) = kin_vis(i,k)*10^(-6);

    elseif p(i,k) >= 50 && p(i,k) <= 1000
        
        % X: Temperature in °C , Y: pressure in bar
        X1= [-150 -125	-100	-75	-50	-25	0	25	50	75	100	125	150	175	200	250	300	400	500	600	700	800	900	1000];
        Y1= [1;5;10;20;30;40;50;60;70;80;90;100;150;200;250;300;350;400;450;500;600;700;800;900;1000];

        X_Kord=[X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;X1;];                     %X-Values must be an array with the size of z
        Y_Kord= [Y1	Y1	Y1	Y1	Y1	Y1	Y1	Y1	Y1	Y1	Y1	Y1	Y1	Y1	Y1	Y1 Y1	Y1	Y1	Y1	Y1	Y1	Y1	Y1];  %Y-Values must be an array with the size of z

        % Z:kinematic viscosity in 10^-7 m/s(values from VDI-Wärmeatlas (messured values)D2.2 Tabelle 14)
        Z_Kord= [30.29	43.37	58.34	75.08	93.49	113.5	135.0	157.9	182.2	207.8	234.6	262.7	292.0	322.4	353.9	420.3	490.7	643.5	811.2	993.0	1188.3	1396.7	1617.8	1851.4
        5.830	8.517	11.56	14.95	18.67	22.70	27.03	31.63	36.52	41.66	47.05	52.69	58.56	64.66	70.99	84.28	98.41	129.0	162.6	199.0	238.1	279.8	324.0	370.8
        2.769	4.165	5.722	7.444	9.326	11.36	13.54	15.86	18.32	20.90	23.61	26.44	29.39	32.45	35.62	42.29	49.37	64.70	81.51	99.74	119.3	140.2	162.3	185.7
        1.217	1.997	2.813	3.701	4.663	5.699	6.806	7.981	9.223	10.53	11.89	13.32	14.81	16.35	17.94	21.29	24.85	32.55	40.98	50.11	59.92	70.38	81.47	93.17
        0.7576	1.281	1.854	2.463	3.118	3.821	4.569	5.362	6.198	7.076	7.995	8.953	9.950	10.98	12.05	14.30	16.68	21.83	27.47	33.57	40.12	47.11	54.51	62.33
        0.7873	0.9288	1.383	1.852	2.353	2.887	3.456	4.057	4.690	5.354	6.049	6.772	7.525	8.305	9.111	10.80	12.60	16.47	20.72	25.31	30.23	35.48	41.04	46.91
        0.8122	0.7293	1.108	1.492	1.899	2.332	2.792	3.278	3.789	4.324	4.884	5.467	6.072	6.700	7.348	8.709	10.15	13.26	16.66	20.34	24.29	28.50	32.95	37.66
        0.8344	0.6278	0.9344	1.258	1.601	1.966	2.353	2.761	3.190	3.640	4.110	4.598	5.106	5.631	6.175	7.314	8.519	11.12	13.96	17.04	20.33	23.84	27.56	31.49
        0.8545	0.6101	0.8198	1.097	1.393	1.708	2.042	2.395	2.766	3.153	3.558	3.980	4.417	4.870	5.338	6.318	7.355	9.593	12.04	14.68	17.51	20.52	23.71	27.08
        0.8733	0.6263	0.7442	0.9813	1.241	1.518	1.812	2.122	2.449	2.790	3.146	3.517	3.902	4.300	4.711	5.572	6.483	8.447	10.59	12.91	15.39	18.03	20.83	23.78
        0.8909	0.6474	0.6960	0.8967	1.126	1.372	1.635	1.912	2.204	2.509	2.828	3.159	3.502	3.857	4.225	4.993	5.806	7.557	9.468	11.53	13.74	16.09	18.58	21.21
        0.9076	0.6677	0.6675	0.8340	1.037	1.259	1.496	1.746	2.010	2.286	2.574	2.873	3.183	3.505	3.836	4.531	5.265	6.845	8.570	10.43	12.42	14.54	16.79	19.16
        0.9821	0.7502	0.6569	0.6968	0.8053	0.9430	1.098	1.2644	1.442	1.628	1.823	2.026	2.236	2.454	2.679	3.150	3.647	4.715	5.878	7.132	8.472	9.896	11.40	12.99
        1.047	0.8143	0.6990	0.6803	0.7297	0.8161	0.9233	1.0441	1.175	1.314	1.461	1.614	1.774	1.939	2.110	2.467	2.844	3.655	4.537	5.486	6.499	7.576	8.714	9.913
        1.107	0.8692	0.7436	0.6971	0.7103	0.7619	0.8370	0.9273	1.028	1.138	1.254	1.377	1.504	1.637	1.775	2.063	2.368	3.023	3.735	4.501	5.318	6.186	7.102	8.067
        1.163	0.9186	0.7857	0.7234	0.7139	0.7413	0.7931	0.8611	0.9407	1.029	1.124	1.225	1.331	1.442	1.557	1.798	2.054	2.605	3.203	3.846	4.532	5.260	6.029	6.837
        1.217	0.9642	0.8251	0.7524	0.7276	0.7376	0.7719	0.8231	0.8862	0.9584	1.038	1.123	1.213	1.307	1.406	1.613	1.833	2.308	2.825	3.380	3.972	4.600	5.263	5.960
        1.268	1.007	0.8623	0.7819	0.7461	0.7432	0.7640	0.8018	0.8520	0.9114	0.9781	1.051	1.128	1.210	1.296	1.477	1.670	2.088	2.543	3.032	3.554	4.106	4.689	5.302
        1.319	1.048	0.8977	0.8109	0.7668	0.7539	0.7641	0.7911	0.8307	0.8798	0.9363	0.9988	1.066	1.138	1.213	1.374	1.546	1.918	2.325	2.762	3.229	3.723	4.244	4.792
        1.368	1.087	0.9315	0.8393	0.7885	0.7677	0.7693	0.7875	0.8183	0.8588	0.9068	0.9608	1.020	1.083	1.150	1.294	1.448	1.784	2.152	2.547	2.969	3.417	3.888	4.384
        1.464	1.162	0.9955	0.8940	0.8328	0.8001	0.7885	0.7931	0.8102	0.8370	0.8715	0.9122	0.9581	1.008	1.062	1.179	1.306	1.586	1.895	2.227	2.583	2.959	3.356	3.773
        1.558	1.234	1.056	0.9461	0.8768	0.8352	0.8140	0.8085	0.8153	0.8318	0.8561	0.8866	0.9224	0.9625	1.006	1.103	1.211	1.450	1.714	2.002	2.308	2.634	2.977	3.337
        1.650	1.303	1.113	0.9959	0.9198	0.8711	0.8424	0.8291	0.8279	0.8363	0.8524	0.8749	0.9027	0.9349	0.9709	1.052	1.144	1.351	1.582	1.834	2.105	2.392	2.694	3.012
        1.742	1.371	1.169	1.044	0.9616	0.9070	0.8723	0.8526	0.8449	0.8466	0.8562	0.8721	0.8933	0.9190	0.9485	1.017	1.096	1.277	1.482	1.707	1.948	2.204	2.475	2.760
        1.833	1.437	1.223	1.090	1.002	0.9426	0.9026	0.8777	0.8645	0.8608	0.8648	0.8752	0.8909	0.9111	0.9352	0.9930	1.061	1.221	1.404	1.606	1.824	2.056	2.301	2.559];


        F = griddedInterpolant(X_Kord', Y_Kord', Z_Kord','spline'); % interpolation 

        nu(i,k) = F({T(i,k),p(i,k)})*10^(-7);     % gives back the interpolated nu-value
        
    else
        nu(i,k) = 1*10^(-7);     % without a nu-value it's difficult to find a solution -> default value
        % errordlg('To calculate the viscosity pressure must be between 1 bar and 1000 bar')    
    end
    end    
end
end

