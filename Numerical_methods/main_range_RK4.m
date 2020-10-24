% Aircraft Design - Long Range Private Jet
% Weight variation with distance

clc;
clear all;
close all;

%% 1. Definition of Constants, Parameters and Variables

% 1.1. CONSTANTS   
Ru = 8.31432;           % Universal Constant for Ideal Gases    [J/mole*K]
% Earth
g = 9.80665;            % Acceleration at Earth's surface       [m/s^2]
T0 = 288.15;            % US Standard Sea Level Temperature     [K]
P0 = 101325;            % Pressure at Sea Level                 [Pa]
Mm = 28.9644*10^-3;     % Molecular Mass                        [kg*mole^-1]    
H_layer = 1e3*[0 11 20 32 47 52 61 69 79 90 100 110 117.776];   % Earth's atmospheric layers
lambda = 1e-3*[-6.5 0 1 2.8 0 -2 -4 -3 0 2 4.36 16.4596 0];     % Earth's atmospheric layers altitude thermal gradient [k/m]
gamma = 1.4;            % Earth's air specific heats relation   [adim]
R = Ru/Mm;              % Gas constant for Earth's air

% 1.2 AIRCRAFT VALUES
c_t_specific = 18/(10^6);       % Kg/(N·s)
c_t_engine = c_t_specific *g;   % 1/s
c_t = c_t_engine*2;             % 2 engines
M = 0.8;                        % Mach Number
C_D0 = 0.016;                   % Parasitic drag constant
k = 0.0045;                     % Lift induced drag constant
S = 50;                         % Wing surface      [m^2]
H = 11000;                      % Cruise altitude   [m]
W_0 = 30000;                    % MTOW              [kg]
Range = (7000 * 1.852)*10^3;    % Aircraft range    [m]

%% 2. PHYSICAL DATA

% Compute base temperatures and pressures
[Tb, Pb] = getBaseTemperaturePressure(R, g, T0, P0, Mm, H_layer, lambda);
% Compute temperature and pressure at cruise altitude
T_h = getTemperatureV2(Tb, H_layer, lambda, H);
P_h = getPressureV2(Tb, Pb, H_layer, lambda, R, g, Mm, H);
rho = getDensityV2(Tb, Pb, H_layer, lambda, R, g, Mm, H);

%% 3. RUNGE-KUTTA RK4

% 3.1 Numerical data
DeltaX = 100;

% Declaration of solution vectors
X_sol = 0:DeltaX:Range;
W_sol = zeros(1, length(X_sol));

% 3.2 Initial conditions
W_sol(1) = W_0;

% 3.3 Runge-Kutta RK4
for i = 1:length(X_sol)-1
    % Functions
    F1 = @(X, W) -c_t/(M*sqrt(gamma*R*T_h))*(0.5*C_D0*rho*(M^2*gamma*R*T_h)*S+k*(2*W^2)/(rho*(M^2*gamma*R*T_h)*S));
    % Computation of coefficients sub 1
    i1 = F1(X_sol(i), W_sol(i));
    % Computation of coefficients sub 2
    i2 = F1(X_sol(i)+ DeltaX/2, W_sol(i)+i1*DeltaX/2);
    % Computation of coefficients sub 3
    i3 = F1(X_sol(i)+ DeltaX/2, W_sol(i)+i2*DeltaX/2);
    % Computation of coefficients sub 4
    i4 = F1(X_sol(i)+ DeltaX, W_sol(i)+i3*DeltaX);
    % Compute next step
    W_sol(i+1) = W_sol(i) + (DeltaX/6)*(i1 + 2*i2 + 2*i3 + i4);
end

%% 4. PLOT SOLUTION

figure(1);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

plot(X_sol/1000,W_sol)
ylim([0,W_0])
xlabel("Range (km)")
ylabel("Weight (kg)")
grid minor;
box on;

% load('datos.mat')
% hold on;
% plot(X_sol_RK/1000,W_sol_RK);

% Fraction between final and initial mass
frac = W_sol(end)/W_0;



