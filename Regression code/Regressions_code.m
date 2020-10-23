% Code for obtaining auxiliar parameters involved in the weight
% determination

% 2020, Aircraft Design

% Authors: 
% Cristian Asensio García
% Juan Garrido Moreno
% Yi Qiang Ji Zhang
% Alexis Leon Delgado
% Alba Molina Cuadrado
% David Morante Torra
% Teresa Peña Mercadé
% Ferran Rubio Vallhonrat
% Iván Sermanoukian Molina
% Santiago Villarroya Calavia

% PREAMBLE 
clear
clc
close all
format long

%% DATA INPUT

% Fuselage dimensions Excel reading
fuselage_length=[xlsread('JET','Sheet1','W6');xlsread('JET','Sheet1','W8:W11');xlsread('JET','Sheet1','W13');xlsread('JET','Sheet1','W20');xlsread('JET','Sheet1','W22')];
fuselage_width=[xlsread('JET','Sheet1','T6');xlsread('JET','Sheet1','T8:T11');xlsread('JET','Sheet1','T13');xlsread('JET','Sheet1','T20');xlsread('JET','Sheet1','T22')];
fuselage_height=[xlsread('JET','Sheet1','U6');xlsread('JET','Sheet1','U8:U11');xlsread('JET','Sheet1','U13');xlsread('JET','Sheet1','U20');xlsread('JET','Sheet1','U22')];

% Operative masses Excel reading
MTOW=[xlsread('JET','Sheet1','N6');xlsread('JET','Sheet1','N8:N11');xlsread('JET','Sheet1','N13');xlsread('JET','Sheet1','N20');xlsread('JET','Sheet1','N22')];
OEW=[xlsread('JET','Sheet1','P6');xlsread('JET','Sheet1','P8:P11');xlsread('JET','Sheet1','P13');xlsread('JET','Sheet1','P20');xlsread('JET','Sheet1','P22')];
engine_mass=[xlsread('JET','Sheet1','X6');xlsread('JET','Sheet1','X8:X11');xlsread('JET','Sheet1','X13');xlsread('JET','Sheet1','X20');xlsread('JET','Sheet1','X22')];

% Baggage volume Excel reading
baggage_volume=[xlsread('JET','Sheet1','S6');xlsread('JET','Sheet1','S8:S11');xlsread('JET','Sheet1','S13');xlsread('JET','Sheet1','S20');xlsread('JET','Sheet1','S22')];

% Estimation of our aircraft fuselage data
b_f_est=mean(fuselage_width);
l_f_est=mean(fuselage_length);

%% LINEAR REGRESSION OF THE BAGGAGE COMPARTMENT VOLUME

% Computation of the product of the squared fuselage width times the length
baggage_volume_abscissa=(fuselage_width.^2).*fuselage_length;

% Regression of a linear function that crosses the coordinate origin
kb_regression=fitlm(baggage_volume_abscissa,baggage_volume,'Intercept',false);

% Graphic plotting 
fig1=figure(1);
set(fig1,'Renderer', 'painters', 'Position', [400 400 500 350]);
hold on

% Scatter plot
scatter(baggage_volume_abscissa,baggage_volume,'d','r')

% Regression line plot
plot(linspace(50,140,20),linspace(50,140,20)*kb_regression.Coefficients{1,1},'b')

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',10)
ylabel('$V_b\,\left[\mathrm{m}^3\right]$','interpreter','latex','FontSize',12)
xlabel('$b_f^2 l_f \left[\mathrm{m}^3\right]$','interpreter','latex','FontSize',12)
 
% Grid format
grid on
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;


%% LINEAR REGRESSION OF THE DELTA W_E
 
% Computation of the correspondent abscissa
Delta_W_e_abscissa=fuselage_length.*(fuselage_width+fuselage_height)/2;

% Delta W_e computation
Delta_W_e=OEW-0.2*MTOW-500-2*engine_mass;
% In the tri-engine case, another engine mass must be added
Delta_W_e(2)=Delta_W_e(2)+engine_mass(2);
 
% Regression of a standard linear function
%Delta_W_e_regression=polyfit(Delta_W_e_abscissa,Delta_W_e,1);
Delta_W_e_regression=polyfit(log10(Delta_W_e_abscissa),log10(Delta_W_e),1);
 
% Graphic plotting 
fig2=figure(2);
set(fig2,'Renderer', 'painters', 'Position', [400 400 500 350]);
hold on

% Scatter plot
scatter(Delta_W_e_abscissa,Delta_W_e,'d','r')

% Regression line plot 
plot(linspace(32,52,20),linspace(32,52,20).^Delta_W_e_regression(1)*10^Delta_W_e_regression(2),'b')
 
% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',10)
xlabel('$l_f\frac{b_f +h_f}{2} \,\left[\mathrm{m}^3\right]$','interpreter','latex','FontSize',12)
ylabel('$\Delta W_e\, \left[\mathrm{kg}\right]$','interpreter','latex','FontSize',12)
set(gca,'yscale','log')
set(gca,'xscale','log')
xlim([32 52])
ylim([3500 9000])

% Grid format
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;

%% LINEAR REGRESSION OF THE ALPHA CONSTANT FOR THE SIMILARITIES CRITERIA

% Regression of a linear function that crosses the coordinate origin
alpha_regression=fitlm(MTOW,OEW,'Intercept',false);

% Graphic plotting 
fig3=figure(3);
set(fig3,'Renderer', 'painters', 'Position', [400 400 500 350]);
hold on

% Scatter plot
scatter(MTOW,OEW,'d','r');

% Regression line plot
plot(linspace(10000,35000,20),linspace(10000,35000,20)*alpha_regression.Coefficients{1,1},'b');

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',10)
xlabel('$\mathrm{MTOW}\,\left[\mathrm{kg}\right]$','interpreter','latex','FontSize',12)
ylabel('$\mathrm{OEW}\,\left[\mathrm{kg}\right]$','interpreter','latex','FontSize',12)
 
% Grid format
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;