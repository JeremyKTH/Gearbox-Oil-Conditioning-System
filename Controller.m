clc
%% ARMAX Model at +-40%
num = [5.154e-05, -4.27e-05];  % B
den = [1, -1.989, 0.9893];     % A
Ts = 0.2;
sys = tf(num, den, Ts)

%% PIDN Values
P = 287.9;
I = 11.42;
D = 1346;
N = 1.864;
b = 0.9272;
c = 0.1973;

%% Controllers
z = tf('z', Ts);

% FeedForward
Gff = b*P + I*(Ts)/(z-1) + c*D*((z-1)*N/(z-1+N*Ts));
% FeedBack
Gc = P + I*(Ts)/(z-1) + D*((z-1)*N/(z-1+N*Ts));
% Plant
Gp = (1.7010e-04)*(z+0.7613)/((z-0.9978)*(z+0.7486));

% Full-state Feedback
Gyr = Gff*Gp/(1 + Gc*Gp);
pzmap(Gyr)
Gyr = minreal(Gyr, 1e-2)

%% Step Response
step(Gyr)
disp('Gp')
pole(Gp)
disp('Gyr')
pole(Gyr)

%% Run for Simulink
[T, R] = tfdata(Gff);
[S, R] = tfdata(Gc);
[B, A] = tfdata(Gp);