%__________________________________________________________________________
% Asparuh Stoyanov, Ramin Sedighi
% Projekt Entwicklung und Simulation
% Schwingungen einer Industriewaschmaschine
%__________________________________________________________________________
% Projekt:  Schwingungsanalyse eines vereinfachten Modells einer 
%           industriellen Waschmaschine mit Matlab und Simulink
%___________________________________________________________________
%  Datenfile
%  File: Unwuchterregung_Anlauf.m;  Bestimmung der  Schwingungsamplitude bei Zentrifuge 
%                                   2-Massenmodel
%________________________________________________________________________
	clear,  close all; clc
%__________________________________________________________________________

s = tf('s');
%--------------------------Parameter des Systems--------------------------%
m = 640;
wasser = 25;
mu = 24+wasser;
k1 = 1467055;  
d = 1.404117761184937e+04;
dT = 0.00001;
ru = 0.330;
alpha = 0.1;                                % Winkelbeschleunigung
%----------------------Eigene Frequenz ohne Daempfung---------------------%
omega_0 = sqrt(k1/(m+mu));
f_0 = omega_0/(2*pi);
%------------------------Kreisfrequenz der Anregung-----------------------%
%omega = 0.5*omega_0;
%omega = 5*omega_0;/
n = 1100;                                   %Zentrifuge U/min
omegaN = 2*pi*(n/60);                       % Kreisfrequenz bei Zentrifuge
omega = omegaN;                             %wir wollen Tilgung bei n=1100 U/min
f = omega/(2*pi);
%---------------Parameter des Tilgers (angepasst an omega)----------------%
kT = 6e+06;                                 %Federsteifigkeit Tilger         
mT = kT/(omegaN)^2;                         %Tilgermasse
%mT = 0.9*kT/(omega)^2,                     %Nicht perfekt angepasste Tilgermasse

%----------------------Frequenzbereich der Anregung-----------------------%
fmin = 0;                                   %Anfangsfrequenz des Anlaufs       
fmax = f;                                   % Stationäre Frequenz nach dem Anlauf
%fmax = 2.5*f;                              % Stationäre Frequenz nach dem Anlauf
%---------------------------Matrizen des Systems--------------------------%
Ai = [(m+mu)*s^2+(d+dT)*s+(k1+kT), -dT*s-kT; -dT*s-kT, mT*s^2+dT*s+kT];
Bi = [1, 0]';
HT = Ai\Bi;

[b1T,a1T] = tfdata(HT(1));
b1T = b1T{:};                               % Koeffizienten der Uebertragungsfunktion 
a1T = a1T{:};                               % von Fe bis x

[b2T,a2T] = tfdata(HT(2));
b2T = b2T{:};                               % Koeffizienten der Uebertragungsfunktion 
a2T = a2T{:};                               % von Fe bis xT

f = linspace(0,f_0*5,1000);
Hx = freqs(b1T, a1T, 2*pi*f);
figure(1);    
clf;
subplot(211), plot(f, 20*log10(abs(Hx)));
title('Amplitudengang des Systems mit Tilger von Fe bis x');
xlabel('Hz');    
ylabel('dB');     
grid on;
axis tight;
subplot(212)
plot(f, angle(Hx)*180/pi);
title('Amplitudengang des Systems mit Tilger von Fe bis x');
xlabel('Hz');    
ylabel('Grad');    
grid on;
axis tight;

%------------------------Aufruf der Simulation----------------------------%
dt = 0.2;     Tfinal = 40;
t = 0:dt:Tfinal;
T1 = 20;
alpha = 2*pi*(fmax-fmin)/T1;
f = ((fmin + (fmax-fmin)*t/T1).*(t<=T1) + fmax*(t>T1));

my_options = simset('MaxStep', dt);
sim('Unwucht4_Frequenzgang.slx', [t], my_options);

t = x.time;
x  = x.signals.values;
xT = xT.signals.values;

figure(2);    clf;
subplot(311), plot(t, f);
title(['Frequenz der Anregung; Eigenfrequenz =',...
    num2str(f_0),' Hz']);
xlabel('Zeit in s');    
grid on; 
ylabel('Hz');
La = axis;
hold on;
plot(La(1:2), [f_0, f_0],'--');
hold off;

subplot(312), plot(t, x);
title('Schwingung der Lage der Hauptmasse');
xlabel('Zeit in s');
grid on;

subplot(313), plot(t, xT);
title('Schwingung der Lage der Tilgersmasse');
xlabel('Zeit in s');
grid on;
