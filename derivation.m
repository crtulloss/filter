% Caleb Rees Tulloss
% Jagannaath Shiva Letchumanan
% ELEN E4215 Filters
% Project: Channel Select Bandpass Filter for Transceiver

clear;
close;

%% Specs
omega_center = 10e6 * 2*pi;
B = 2e6 * 2*pi;

omega_p1 = omega_center - B/2;
omega_p2 = omega_center + B/2;

omega_0 = sqrt(omega_p1 * omega_p2);
Q = omega_0 / B;

omega_s1 = 7.5e6 * 2*pi;
omega_s2 = 12.5e6 * 2*pi;

gain_dB = 10;
gain = 10^(gain_dB / 20);

alpha_max = 0.5;   % dB
alpha_min = 60;    % dB
epsilon = sqrt(10^(alpha_max / 10) - 1);

%% Transform to Lowpass for Design

% normalize frequency axis
omega_s1_norm = omega_s1 / omega_0;
omega_s2_norm = omega_s2 / omega_0;
omega_p1_norm = omega_p1 / omega_0;
omega_p2_norm = omega_p2 / omega_0;

Omega_s1 = Q*(omega_s1_norm^2 - 1)/omega_s1_norm;
Omega_s2 = Q*(omega_s2_norm^2 - 1)/omega_s2_norm;

Omega_p1 = Q*(omega_p1_norm^2 - 1)/omega_p1_norm;
Omega_p2 = Q*(omega_p2_norm^2 - 1)/omega_p2_norm;

% Omega_s2_norm is a more strict requirement (sharper transition),
% so use that

atten_for_graph = alpha_min + 20*log(1/epsilon);
% can look at the graph - looks like 5

% also check here
[n,Wp] = ellipord(Omega_p2, Omega_s2, alpha_max, alpha_min, 's');
% confirmation: order 5

%% Lowpass Transfer Function

% look in design table - doesn't have exactly what we want...??

% or use Matlab design - can we just use this?
[ze,pe,ke] = ellip(n,0.5,60,Omega_p2,'s');
[be,ae] = zp2tf(ze,pe,ke);
[he,we] = freqs(be,ae,4096);
figure
plot(we,mag2db(abs(he)));
title('Normalized Lowpass Response: 5^{th}-order Elliptic');
ylabel('|H(j\omega)| (dB)');
xlabel('Normalized \Omega (rad/s)');
