Ts = 0.4;
z = tf('z', Ts);
%% Best +/-40% Model: ARMAX
b0 = 3.264e-05;
a1 = -1.98;
a0 = 0.9805;
B = [b0];         % B
A = [1, a1, a0];  % A
Ts = 0.4;
Gp_armax = tf(B, A, Ts);

%% Controller Values
% Tr = 132.8, Ts = 219.6
P = 0.041945;
I = 0.20972;
D = 0;
N = 0;
b = 1;
c = 1;

% Tr = 12, Ts = 32
% P = 111.4764;
% I = 2.5925;
% D = 936.5568;
% N = 0.62227;
% b = 0.93997;
% c = 0.36959;

Gc_armax = P + I*Ts/(z-1) + D*(N)/(1+(N*Ts)/(z-1));
Gff_armax = b*P + I*Ts/(z-1) + c*D*(N)/(1+(N*Ts)/(z-1));
Gyr_armax = Gff_armax*Gp_armax/(1+Gc_armax*Gp_armax);

figure(1)
stepplot(Gp_armax, Gyr_armax)
legend('Original', 'Controller')
grid on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TF Model
b0 = 0.001273;
a = -0.9965;
B = [b0];    % B
A = [1, a];  % A
Ts = 0.4;

Gp_tf = tf(B, A, Ts);
%% Controller Values
% Tr = 128.4, Ts = 220.8
P = 13.05;
I = 0.2043;
b = 0.01175;

% Tr = 16.04, Ts = 52
% P = 42.8951;
% I = 4.0336;
% b = 0.049622;

Gc_tf = P + I*Ts/(z-1);
Gff_tf = b*P + I*Ts/(z-1);
Gyr_tf = Gff_tf*Gp_tf/(1+Gc_tf*Gp_tf);

figure(2)
stepplot(Gp_tf, Gyr_tf)
legend('Original', 'Controller')
grid on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SS
A_ss = [[1.037, 0.4309];
       [-0.01128, 0.8728]];
B_ss = [-0.003958;
        0.0011];
C_ss = [-3.142, -0.6329];
D_ss = 0;

[B, A] = ss2tf(A_ss, B_ss, C_ss, D_ss);
b1 = B(1); b0 = B(2);
a1 = A(2); a0 = A(3);
Gp_ss = tf(B, A, Ts);

%% Controller Values
% Tr = 123, Ts = 332
% P = 4.6682;
% I = 0.1104;
% b = -128.2122;
% N = 27.7958;
% b = 0.012615;
% c = 0.007922;

% Tr = 15.6, Ts = 94.8
P = 6.4143;
I = 1.2514;
b = -3.1995;
N = 0.71906;
b = 0.042499;
c = 0.0036211;

Gc_ss = P + I*Ts/(z-1) + D*(N)/(1+(N*Ts)/(z-1));
Gff_ss = b*P + I*Ts/(z-1) + c*D*(N)/(1+(N*Ts)/(z-1));
Gyr_ss = Gff_ss*Gp_ss/(1+Gc_ss*Gp_ss);

figure(3)
stepplot(Gp_ss, Gyr_ss)
legend('Original', 'Controller')
grid on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison
figure(4)
stepplot(Gyr_armax, Gyr_tf, Gyr_ss)
legend('ARMAX', 'TF', 'SS')
grid on