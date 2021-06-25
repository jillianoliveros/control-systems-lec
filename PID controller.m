% ECE21122 - Control Systems Theory
% Machine Problem 1
% Jessica Ricci A. Lapeña
% Jillian Eunice B. Oliveros
% Joshua Alexander B. Villa
% 3ECE-C

close all
clear all
clc
syms s zeta zc


% step 1: evaluating the performance of the existing system

% for the uncompensated system
OS = exp(-(pi*zeta)/(sqrt(1-zeta^2)))*100 == 16.3024;
solved_zeta = round(vpasolve(OS,zeta),5)
% the damping ratio is 0.500
num_uncomp = [(s+3.4)];
den_uncomp = [(s+2.5)*(s+4.5)*(s+8)];
num_uncomp = sym2poly(num_uncomp);
den_uncomp =sym2poly(den_uncomp);
Gs_uncomp = tf(num_uncomp,den_uncomp);
% plot the rlocus of the uncompensated system
figure
rlocus(Gs_uncomp)
axis([-10 5 -5 15])
% move the point on the rlocus plot to find the %OS of 16.3%
% the gain (K) is 107
% the dominant poles are at -5.82 + 10.1i
% the damping ratio (zeta) is 0.5
% the oscillation frequency (wn) is 11.6
%substituting K and computing for the closed-loop transfer function
num_uncomp1 = [107*(s+3.4)];
den_uncomp1 = [(s+2.5)*(s+4.5)*(s+8)];
num_uncomp1 = sym2poly(num_uncomp1);
den_uncomp1 =sym2poly(den_uncomp1);
Gs_uncomp1 = tf(num_uncomp1,den_uncomp1)
Ts_uncomp = feedback(Gs_uncomp1,1)
[num, den] = tfdata(Ts_uncomp);
den_fb = cell2mat(den);
roots_den = roots(den_fb)
% the third pole is at -3.3575 based on the roots
% plot the rlocus of the uncompensated system with K
figure
rlocus(Ts_uncomp)
% second-order approximation is valid
% the third pole at -3.3575 and zero at -3.4 will cancel each other
step_info_uncomp = stepinfo(Ts_uncomp)


% step 2: PD controller design

% the dominant poles are at -5.82 +- 10.1i
% peak time within 50% to 75% that of the uncompensated system
% peak time (Tp) = pi/wn
% multiply the imaginary and real parts to 2
% real part of the compensated system = -11.64
% imaginary part of the compensated system = 20.2
% summation of zeros - summation of poles (open-loop transfer function)
s1 = -11.64+20.2i;
% angle contribution = (s+3.4)-[(s+2.5)+(s+4.5)+(s+8)]
% substitute s1 to the equation
z1 = (s1+3.4);
z1_deg = 180-((atan(abs(imag(z1)/real(z1))))*180/pi);
p1 = (s1+2.5);
p1_deg = 180-((atan(abs(imag(p1)/real(p1))))*180/pi);
p2 = (s1+4.5);
p2_deg = 180-((atan(abs(imag(p2)/real(p2))))*180/pi);
p3 = (s1+8);
p3_deg = 180-((atan(abs(imag(p3)/real(p3))))*180/pi);
ang_sum = z1_deg - (p1_deg+p2_deg+p3_deg)
% the sum of the angles to the desired pole is −211.8357°
ang = -(ang_sum) - 180;
ang1 = tand(ang) == (20.2/(zc-11.64));
comp_zero = round(solve(ang1,zc),4)
% the compensator zero is -zc = -44.174
% PD controller transfer function: Gs_PD = s+comp_zero;
num_PD = [(s+comp_zero)*(s+3.4)];
den_PD = [(s+2.5)*(s+4.5)*(s+8)];
num_PD = sym2poly(num_PD);
den_PD =sym2poly(den_PD);
Gs_PD = tf(num_PD,den_PD);
% plot the rlocus of the compensated system
figure
rlocus(Gs_PD)
axis([-120 10 -40 40])
% move the point on the rlocus plot to find the %OS of 16.3%
% the gain (K) is 9.26
% the dominant poles are at -10.4 + 18.1i
% the damping ratio (zeta) is 0.5
% the %OS is at 16.3%
% the oscillation frequency (wn) is 20.9
% substituting K and computing for the closed-loop transfer function
num_PD_comp = [9.26*(s+comp_zero)*(s+3.4)];
den_PD_comp = [(s+2.5)*(s+4.5)*(s+8)];
num_PD_comp = sym2poly(num_PD_comp);
den_PD_comp =sym2poly(den_PD_comp);
Gs_PD_comp = tf(num_PD_comp,den_PD_comp)
Ts_PD = feedback(Gs_PD_comp,1)
[num, den] = tfdata(Ts_PD);
den_fb = cell2mat(den);
roots_den = roots(den_fb)
% the third pole is at  -3.3879 based on the roots
% plot the rlocus of the compensated system with K
figure
rlocus(Ts_PD)
% second-order approximation is valid
% the third pole at  -3.3879 and zero at -3.4 will cancel each other
step_info_PD = stepinfo(Ts_PD)


% step 3: simulation and verification

figure
step(Ts_PD,'b',Ts_uncomp,'k')
legend('PD-compensated','Uncompensated')
% the peak time was reduced from 0.3110s to 0.1736s
% the steady-state error was reduced from 0.1983 to 0.0607


% step 4: PI controller design

% choose a value of zero close to the origin = -0.1
% PD controller transfer function: Gs_PI = (s+0.1)/s
num_PI = [(s+0.1)*(s+comp_zero)*(s+3.4)];
den_PI = [s*(s+2.5)*(s+4.5)*(s+8)];
num_PI = sym2poly(num_PI);
den_PI =sym2poly(den_PI);
Gs_PI = tf(num_PI,den_PI);
% plot the rlocus of the compensated system
figure
rlocus(Gs_PI)
% move the point on the rlocus plot to find the %OS of 16.3%
% the gain (K) is 8.23
% the dominant poles are at -9.88 + 17.1i
% the damping ratio (zeta) is 0.5
% the %OS is at 16.3%
% the oscillation frequency (wn) is 19.8
% substituting K and computing for the closed-loop transfer function
num_PI_comp = [8.23*(s+0.1)*(s+comp_zero)*(s+3.4)];
den_PI_comp = [s*(s+2.5)*(s+4.5)*(s+8)];
num_PI_comp = sym2poly(num_PI_comp);
den_PI_comp =sym2poly(den_PI_comp);
Gs_PI_comp = tf(num_PI_comp,den_PI_comp)
Ts_PI = feedback(Gs_PI_comp,1)
[num, den] = tfdata(Ts_PI);
den_fb = cell2mat(den);
roots_den = roots(den_fb)
% plot the rlocus of the compensated system with K
figure
rlocus(Ts_PI)
% second-order approximation is invalid
% the third pole at -3.3860
% the fourth pole is at -0.0935
step_info_PI = stepinfo(Ts_PI)


% step 5: PID controller design

num_PID_comp = [8.23*(s+0.1)*(s+comp_zero)*(s+3.4)];
den_PID_comp = [s*(s+2.5)*(s+4.5)*(s+8)];
num_PID_comp = sym2poly(num_PID_comp);
den_PID_comp =sym2poly(den_PID_comp);
Gs_PID_comp = tf(num_PI_comp,den_PID_comp)
Ts_PID = feedback(Gs_PI_comp,1)
[num, den] = tfdata(Ts_PID);
den_fb = cell2mat(den);
roots_den = roots(den_fb)
figure
rlocus(Ts_PID)
step_info_PID = stepinfo(Ts_PID)


% step 6: simulation and verification of PID controller 

figure
step(Ts_PID,'r',Ts_PD,'b',Ts_uncomp,'k')
legend('PID-compensated','PD-compensated','Uncompensated')
axis([0 1 0 1.2])
