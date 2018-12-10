% Caleb Rees Tulloss
% Jagannaath Shiva Letchumanan
% ELEN E4215 Filters
% Project: Channel Select Bandpass Filter for Transceiver

clear;
close all;

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

%% Transform to Lowpass for Design

% normalize frequency axis
omega_s1_norm = omega_s1 / omega_0;
omega_s2_norm = omega_s2 / omega_0;
omega_p1_norm = omega_p1 / omega_0;
omega_p2_norm = omega_p2 / omega_0;

% transform to lowpass
Omega_s1 = Q*(omega_s1_norm^2 - 1)/omega_s1_norm;
Omega_s2 = Q*(omega_s2_norm^2 - 1)/omega_s2_norm;
Omega_p1 = Q*(omega_p1_norm^2 - 1)/omega_p1_norm;
Omega_p2 = Q*(omega_p2_norm^2 - 1)/omega_p2_norm;

% Omega_s2_norm is a more strict requirement (sharper transition),
% so use that

epsilon = sqrt(10^(alpha_max / 10) - 1);
atten_for_graph = alpha_min + 20*log(1/epsilon);
% can look at the graph - looks like 5

% also check here
[n,Wp] = ellipord(Omega_p2, Omega_s2, alpha_max, alpha_min, 's');
% confirmation: order 5

%% Lowpass Transfer Function

% or use Matlab design
[ze,pe,ke] = ellip(n,0.5,60,Omega_p2,'s');
% don't forget gain!
ke = ke*gain;

% plot lowpass just for fun
[be,ae] = zp2tf(ze,pe,ke);
[he,we] = freqs(be,ae,4096);
figure
plot(we,mag2db(abs(he)));
title('Normalized Lowpass Response: 5^{th}-order Elliptic');
ylabel('|H(j\Omega)| (dB)');
xlabel('Normalized \Omega (rad/s)');

%% Transform to Bandpass
% note that the textbook derives these for omega_0=1, while
% we use our actual omega_0 to arrive at the final transfer function

% calculate w0 and Q for each of the cc pole pairs
p1_w0 = abs(pe(1));
p1_Q = abs(p1_w0 / (2*real(pe(1))));
p2_w0 = abs(pe(3));
p2_Q = abs(p2_w0 / (2*real(pe(3))));

% follow the textbook derivation p.385-386 to calculate the
% bandpass w0 and Q - 2 cc pairs for each of the lowpass cc pairs
p1_a = p1_w0 / p1_Q;
p1_b = p1_w0^2;
p2_a = p2_w0 / p2_Q;
p2_b = p2_w0^2;

p1_q = sqrt(Q/p1_a * ((2*Q/p1_a + p1_b/(2*p1_a*Q)) + sqrt(...
    (2*Q/p1_a + p1_b/(2*p1_a*Q))^2 - 1)));
p2_q = sqrt(Q/p2_a * ((2*Q/p2_a + p2_b/(2*p2_a*Q)) + sqrt(...
    (2*Q/p2_a + p2_b/(2*p2_a*Q))^2 - 1)));

p1_w01 = omega_0*(p1_a*p1_q/(2*Q)+0.5*sqrt(p1_b/(Q^2) - 1/(p1_q^2)));
p1_w02 = omega_0^2/p1_w01;

p2_w01 = omega_0*(p2_a*p2_q/(2*Q)+0.5*sqrt(p2_b/(Q^2) - 1/(p2_q^2)));
p2_w02 = omega_0^2/p2_w01;

% for the lowpass real pole, the bandpass cc pole pair
% has w0=omega_0. calculate Q
p3 = abs(pe(5));
p3_w0 = omega_0;
p3_q = Q/p3;

% follow the textbook derivation p. 387-388 to calculate the bandpass zeros
% signs are flipped in 2 and 4 because -a * -1/a = 1 as well as a/a=1!
% ensures that zeros are cc pairs
z1 = imag(ze(1));
z1_w1 = omega_0*(sqrt(z1^2/(4*Q^2) + 1) - z1/(2*Q));
z1_w2 = omega_0*(sqrt(z1^2/(4*Q^2) + 1) + z1/(2*Q));
z2 = imag(ze(2));
z2_w1 = omega_0*(-sqrt(z2^2/(4*Q^2) + 1) + z2/(2*Q));
z2_w2 = omega_0*(-sqrt(z2^2/(4*Q^2) + 1) - z2/(2*Q));
z3 = imag(ze(3));
z3_w1 = omega_0*(sqrt(z3^2/(4*Q^2) + 1) - z3/(2*Q));
z3_w2 = omega_0*(sqrt(z3^2/(4*Q^2) + 1) + z3/(2*Q));
z4 = imag(ze(4));
z4_w1 = omega_0*(-sqrt(z4^2/(4*Q^2) + 1) + z4/(2*Q));
z4_w2 = omega_0*(-sqrt(z4^2/(4*Q^2) + 1) - z4/(2*Q));

% check the transfer function - need to put in pole-zero form
% zero at origin from transformation of 1st-order pole
bandpass_z = [0 z1_w1 z1_w2 z2_w1 z2_w2 z3_w1 z3_w2 z4_w1 z4_w2]';
bandpass_z = 1i*bandpass_z;

% for the purpose of matlab filter plotting, which likes explicit cc poles
% instead of w0,Q for a pole pair
a = 1;
b = p1_w01/p1_q;
c = p1_w01^2;
p = [a b c];
bp12 = roots(p);

a = 1;
b = p1_w02/p1_q;
c = p1_w02^2;
p = [a b c];
bp34 = roots(p);

a = 1;
b = p2_w01/p2_q;
c = p2_w01^2;
p = [a b c];
bp56 = roots(p);

a = 1;
b = p2_w02/p2_q;
c = p2_w02^2;
p = [a b c];
bp78 = roots(p);

a = 1;
b = p3_w0/p3_q;
c = p3_w0^2;
p = [a b c];
bp910 = roots(p);

bandpass_p = vertcat(bp12, bp34, bp56, bp78, bp910);

% factor of omega_0/Q from the 1st order pole
bandpass_k = ke * omega_0 / Q;

% plot normalzied bandpass
[bandpass_b,bandpass_a] = zp2tf(bandpass_z,bandpass_p,bandpass_k);
[bandpass_h,bandpass_w] = freqs(bandpass_b,bandpass_a,10000);
figure
plot(bandpass_w/(1e6*2*pi),mag2db(abs(bandpass_h)));
title('Normalized Bandpass Response: 10^{th}-order Elliptic');
ylabel('|H(j\omega)| (dB)');
xlabel('f (MHz)');

%% Pole/zero pairing
% numx defines the numerator term (s^2 + numx)

% zero at origin with p1_w02
num1 = 0;

% zero at +/- j4.71e7
num2 = bandpass_z(2) * bandpass_z(5);
% with p2_w02

% zero at +/- j5.19e7
num3 = bandpass_z(6) * bandpass_z(9);
% with p3_w0

% zero at +/- j7.52e7
num4 = bandpass_z(7) * bandpass_z(8);
% with p2_w01

% zero at +/- j8.29e7
num5 = bandpass_z(3) * bandpass_z(4);
% with p1_w01

%% TF ordering

%Q factor of first set of poles
Q1 = abs(abs(bp12(1))/(2*real(bp12(1))));

%Q factor of second set of poles
Q2 = abs(abs(bp34(1))/(2*real(bp34(1))));

%Q factor of third set of poles
Q3 = abs(abs(bp56(1))/(2*real(bp56(1))));

%Q factor of fourth set of poles
Q4 = abs(abs(bp78(1))/(2*real(bp78(1))));

%Q factor of fifth set of poles
Q5 = abs(abs(bp910(1))/(2*real(bp910(1))));

%Q5 is the least, followed by Q4 (Q3 has same Q) and Q2 (Q1 has the same Q)
%The order will be p3_w0, p2_w02, p2_w01, p1_w02 and p1_w01.
%Re-ordering the zero and pole vectors
bandpass_p2 = flipud(bandpass_p);
bandpass_z2 = [bandpass_z(6) bandpass_z(9) bandpass_z(2) bandpass_z(5)...
    bandpass_z(7) bandpass_z(8) bandpass_z(1) bandpass_z(1) bandpass_z(3)...
    bandpass_z(4)]';

%% Gain splitting

%Synthesizing individual biquad transfer functions
[t1_b,t1_a] = zp2tf(bandpass_z2(1:2),bandpass_p2(1:2),1);
T1 = tf(t1_b,t1_a);
[t2_b,t2_a] = zp2tf(bandpass_z2(3:4),bandpass_p2(3:4),1);
T2 = tf(t2_b,t2_a);
[t3_b,t3_a] = zp2tf(bandpass_z2(5:6),bandpass_p2(5:6),1);
T3 = tf(t3_b,t3_a);
[t4_b,t4_a] = zp2tf(bandpass_z2(7:8),bandpass_p2(7:8),1);
T4 = tf(t4_b,t4_a);
[t5_b,t5_a] = zp2tf(bandpass_z2(9:10),bandpass_p2(9:10),1);
T5 = tf(t5_b,t5_a);

%Maximum values of the biquads added one at a time
M1 = getPeakGain(T1);
M2 = getPeakGain(T1*T2);
M3 = getPeakGain(T1*T2*T3);
M4 = getPeakGain(T1*T2*T3*T4);
M5 = getPeakGain(T1*T2*T3*T4*T5);

%Optimal gains assigned to each biquad
k(1) = bandpass_k*M5/M1;
k(2) = M1/M2;
k(3) = M2/M3;
k(4) = M3/M4;
k(5) = M4/M5;