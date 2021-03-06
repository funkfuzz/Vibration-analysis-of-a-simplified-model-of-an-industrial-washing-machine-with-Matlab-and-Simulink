%__________________________________________________________________________
% Asparuh Stoyanov, Ramin Sedighi
% Projekt Entwicklung und Simulation
% Schwingungen einer Industriewaschmaschine
%__________________________________________________________________________
% Projekt:  Schwingungsanalyse eines vereinfachten Modells einer 
%           industriellen Waschmaschine mit Matlab und Simulink
%___________________________________________________________________
%  Datenfile
%  File: Unwucterregung_mit_Tilger.m;   Bestimmung der  Schwingungsamplitude 
%                                       bei Zentrifuge 2-Massenmodel                                    
%________________________________________________________________________
	clear,  close all; clc
%__________________________________________________________________________
clear;
s = tf('s');                                    % Variable s für die Übertragungsfunktion
%---------------------------- Systemparameter ----------------------------%
m = 640;                                    %Nettogewicht Waschmaschine
wasser = 25;                                %Wassergewicht
mu = 24+wasser;                             %Masse Unwucht Wäsche+Wasser
k1 = 1467055;                               %Federsteifigkeit Waschmaschine
d = 1.404117761184937e+04;                  %Dämpfung 
ru = 0.330;                                 %Radius zur Unwuchtmassewasser
dT = 0.001;                                 %Dämpfung Tilger
% ----------------------Eigenfrequenz ohne Daempfung----------------------%
omega_0 = sqrt(k1/(m+mu));
f_0 = omega_0/(2*pi);
% -----------------------Kreisfrequenz der Anregung-----------------------%
%omega = 0.5*omega_0;
%omega = 5*omega_0;/
n = 1100; %Zentrifuge U/min
omegaN = 2*pi*(n/60);                       % Kreisfrequenz bei Zentrifuge
omega = omegaN;                             %wir wollen Tilgung bei n=1100 U/min
f = omega/(2*pi);

%---------------Parameter des Tilgers (angepasst an omega)----------------%
kT = 6e+06;                                 %Federsteifigkeit Tilger         
mT = kT/(omegaN)^2;                         %Tilgermasse
%mT = 0.85*kT/(omega)^2,                     %Nicht perfekt angepasste Tilgermasse

%------------------Amplitude Hauptmasse bei omega ohne Tilger-------------%

H1 = 1/((m+mu)*s^2+d*s+k1);                 %Übertragungsfunktion 
[b1,a1] = tfdata(H1);
b1 = b1{:};                                 % Koeffizienten der Übertragungsfunktion
a1 = a1{:};
zaehler1 = polyval(b1, 1i*omega);
nenner1  = polyval(a1, 1i*omega);
ampl_x = mu*ru*omega^2*abs(zaehler1/nenner1); % Amplitude ohne Tilger

%-----------------Übertragungsfunktionen mit Tilger-----------------------%
%Systemmatrizen
Ai = [(m+mu)*s^2+(d+dT)*s+(k1+kT), -dT*s-kT; -dT*s-kT, mT*s^2+dT*s+kT];
Bi = [1, 0]';
HT = Ai\Bi;
%--------------Amplitude Waschmaschine x bei omega mit Tilger-------------%
[b1T,a1T] = tfdata(HT(1));
b1T = b1T{:};                               % Koeffizienten Übertragungsfunktion 
a1T = a1T{:};                               % von Fe zu x
zaehler1T = polyval(b1T, 1i*omega);
nenner1T  = polyval(a1T, 1i*omega);
ampl_xT = mu*ru*omega^2*abs(zaehler1T/nenner1T); % Ampl. mit Tilger

%-------------------Amplitude Tilger xT bei omega-------------------------%
[b2T,a2T] = tfdata(HT(2));
b2T = b2T{:};                              % Koeffizienten Übertragungsfunktion 
a2T = a2T{:};                              % von Fe zu xT
zaehler2T = polyval(b2T, 1i*omega);
nenner2T  = polyval(a2T, 1i*omega);
ampl_xmT = mu*ru*omega^2*abs(zaehler2T/nenner2T); % Amplitude des Tilgers

gewinn = ampl_x/ampl_xT                    % Gewinn

%-------Amplitudenfunktion von x ohne Tilger abhängig von omega-----------%
f_t = linspace(0,f*2,10000);
omega_t = 2*pi*f_t;
zaehler1 = polyval(b1, 1i*omega_t);
nenner1  = polyval(a1, 1i*omega_t);
ampl_x = mu*ru*omega_t.^2.*abs(zaehler1./nenner1);
figure(1);  
subplot(211), plot(f_t, ampl_x);
title('Amplitude x ohne Tilger');
xlabel('Hz');
grid on;
hold on;     
La = axis;
plot([f, f], [La(3), La(4)],'r'); 
hold off;
subplot(212), plot(f_t, ampl_x);
title('Amplitude x ohne Tilger (Ausschnitt)');
xlabel('Hz');
grid on;
hold on;
La = axis;
plot([f, f], [La(3), La(4)],'g');
hold off;
axis([f*0.7, f*1.3, La(3:4)]);

%----------Amplitudenfunktion von x mit Tilger abhängig von omega---------%
zaehler1T = polyval(b1T, 1i*omega_t);
nenner1T  = polyval(a1T, 1i*omega_t);
ampl_xT = mu*ru*omega_t.^2.*abs(zaehler1T./nenner1T);
figure(2);  
subplot(211), plot(f_t, ampl_xT);
title('Amplitude x mit Tilger');
xlabel('Hz');    
grid on;
hold on;    
La = axis;
plot([f, f], [La(3), La(4)],'r');  
hold off;
subplot(212), plot(f_t, ampl_xT);
title('Amplitude x mit Tilger (Ausschnitt)');
xlabel('Hz');     
grid on;
hold on;  
La = axis;
plot([f, f], [La(3), La(4)],'r');   
hold off;
axis([f*0.7, f*1.3, La(3:4)]);

%----------Amplitudenfunktion von xT des Tilgers abhängig von omega-------%
zaehler2T = polyval(b2T, 1i*omega_t);
nenner2T  = polyval(a2T, 1i*omega_t);
ampl_xmT = mu*ru*omega_t.^2.*abs(zaehler2T./nenner2T);
figure(3);  
subplot(211)
plot(f_t, ampl_xmT);
title('Amplitude xT des Tilgers');
xlabel('Hz');    
grid on;
hold on;   
La = axis;
plot([f, f], [La(3), La(4)],'r');  
hold off;
subplot(212)
plot(f_t, ampl_xmT);
title('Amplitude xT des Tilgers (Ausschnitt)');
xlabel('Hz');   
grid on;
hold on;    
La = axis;
plot([f, f], [La(3), La(4)],'r');   
hold off;
axis([f*0.7, f*1.3, La(3:4)]);










