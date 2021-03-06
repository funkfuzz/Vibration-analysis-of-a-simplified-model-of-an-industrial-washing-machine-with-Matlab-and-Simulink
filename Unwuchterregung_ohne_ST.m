%__________________________________________________________________________
% Asparuh Stoyanov, Ramin Sedighi
% Projekt Entwicklung und Simulation
% Schwingungen einer Industriewaschmaschine
%__________________________________________________________________________
% Projekt:  Schwingungsanalyse eines vereinfachten Modells einer 
%           industriellen Waschmaschine mit Matlab und Simulink
%___________________________________________________________________
%  Datenfile
%  File: Unwuchterregung_ohne_ST.m;  Bestimmung der  Schwingungsamplitude bei Zentrifuge 
%                                    1-Massenmodel
%________________________________________________________________________
	clear,  close all; clc
%__________________________________________________________________________

%---------------------------- Systemparameter ----------------------------%
m = 640;                                    %Nettogewicht Waschmaschine
wasser = 25;                                %Wassergewicht
mu = 24+wasser;                             %Masse Unwucht Wäsche+Wasser
k1 = 1467055;                               %Federsteifigkeit Waschmaschine
d = 1.404117761184937e+04;                  %Dämpfung 
ru = 0.330;                                 %Radius zur Unwuchtmasse
% -------- Amplitude x
omega_0 = sqrt(k1/(m+mu));                  % Eigenkreisfrequenz ohne Dämpfung 
f_0 = omega_0/(2*pi);                       % Eigenfrequenz ohne Dämpfung
fmin = f_0/100;                             % minimale Frequenz              
fmax = f_0*50;                              % maximale Frequenz
a1 = round(log10(fmin));                    % Minimalwert für log-Vektor f 
a2 = round(log10(fmax));                    % Maximalwert für log-Vektor f
f = logspace(a1, a2, 1000);                 % Logspace-Vektor für logarithmischer Plot
omega = (2*pi)*f;
ampl_x = mu*(omega.^2)*ru./sqrt((k1-(m+mu)*... % amplitude X
    (omega.^2)).^2 + (d*omega).^2);

%-------------------------------Plots-------------------------------------%
clf;
figure(1);    
semilogx(f, ampl_x);                        % Plots die Amplitude logarithmisch
title('Amplitude x');
xlabel('Hz');    
grid on;      
ylabel('ampl_x');    
hold on;
La = axis;    
loglog([f_0, f_0], [La(3), La(4)], 'g');
hold off;