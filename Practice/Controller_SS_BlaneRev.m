clc

T_train = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\Ts = 0.4\Training_Data\D_48404_train.csv');
T_train.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1'};
% T_train = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\Ts = 1\Training_Data\D_48404_train.csv');
% T_train.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11'};

T_test = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\Ts = 0.4\Testing_Data\D_48404_test.csv');
T_test.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1'};
% T_test = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\Ts = 1\Testing_Data\D_48404_test.csv');
% T_test.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11'};

u_train = table2array(T_train(:, 1));   % OTSV1
y_train = table2array(T_train(:, 2));   % TV12
u_test =  table2array(T_test(:, 1));   % OTSV1
y_test =  table2array(T_test(:, 2));   % TV12

Ts = .4; %.4 For system ID
Ts2 = .4; %1 After System ID
Ts3 = 1; % For Simulink

%% Training Data
data_train = iddata(y_train, u_train, Ts);
data_train.Name = 'Data_Train';
data_train.TimeUnit = 'seconds';
data_train.InputName = 'OTSV1';   data_train.InputUnit = 'Percentage';
data_train.OutputName = 'TV12';   data_train.OutputUnit = 'Celsius';

%% Testing Data
data_test = iddata(y_test, u_test, Ts);
data_test.Name = 'Data_Test';
data_test.TimeUnit = 'seconds';
data_test.InputName = 'OTSV1';   data_test.InputUnit = 'Percentage';
data_test.OutputName = 'TV12';   data_test.OutputUnit = 'Celsius';

%% Detrending Data
T_train = getTrend(data_train);
T_test = getTrend(data_test);
T_train.InputOffset = 40;
T_train.OutputOffset = 48;
T_test.InputOffset = 40;
T_test.OutputOffset = 48;

data_train = detrend(data_train, T_train);
data_test = detrend(data_test, T_test);

cost_func = 'NRMSE';

%% SS System Identification

sysSS = n4sid(data_train, 2, 'Ts', Ts);
%sysSS = ssregest(data_train, 2, 'Ts', Ts);

opt = compareOptions('InitialCondition','e');
figure(1)
compare(data_test, sysSS, opt)

Gp_origin = tf(sysSS);
Gp_cont = d2c(Gp_origin, 'zoh'); %converting to cont.

[num, den] = tfdata(Gp_origin); %2 poles, 1 zero
b2 = num{1}(1); %z^2 (0)
b1 = num{1}(2); %zero z^1
b0 = num{1}(3); %zero z^0

a2 = den{1}(1); %z^2 (1)
a1 = den{1}(2); %pole z^1
a0 = den{1}(3); %pole z^0

B = [b2, b1, b0];  % B
A = [a2, a1, a0];  % A

%Changing sampling time for later
Gp = tf(B, A, Ts2)

poles_Gp_disc = pole(Gp)
zeros_Gp_disc = zero(Gp)

figure(2)
pzmap(Gp)
 
%Gp = c2d(Gp_cont, Ts, 'zoh')

%% Choose Poles (w_m, zeta_m) (w_o, zeta_o)

%////////////////////////// A_m /////////////////////////////

% -----------------FROM CONTINUOUS POLES:--------------------

w_m = .028; % (.028 IS SLOWEST w/ .7 zeta) (.020 IS SLOWEST w/ 1 zeta)
zeta_m = 1; %

% z^2 
p1_m = 2*exp(-zeta_m*w_m*Ts2)*cos(w_m*Ts2*sqrt(1-(zeta_m)^2)); %
p0_m = exp(-2*zeta_m*w_m*Ts2);

syms z1m z2m
eqn1 = 2*exp(-zeta_m*w_m*Ts2)*cos(w_m*Ts2*sqrt(1-(zeta_m)^2)) == z1m + z2m;
eqn2 = exp(-2*zeta_m*w_m*Ts2) == z1m*z2m;

sol = solve([eqn1, eqn2], [z1m, z2m]);
z1m = double(sol.z1m)
z2m = double(sol.z2m)

% -----------------DIRECT DISCRETE POLES:-----------------------

% z1m = .991
% z2m = .781
% 
% p1_m = z1m + z2m; 
% p0_m = z1m*z2m;

%////////////////////////// A_O /////////////////////////////

% ----------------FROM CONTINUOUS POLES:---------------------

w_o = .0235; % (.02 is smallest) ()
zeta_o = 1;

% Den: z^2 - (p1_o)*z + (p0_o) == 0
p1_o = 2*exp(-zeta_o*w_o*Ts2)*cos(w_o*Ts2*sqrt(1-(zeta_o)^2));
p0_o = exp(-2*zeta_o*w_o*Ts2);

syms z1o z2o
eqn1 = 2*exp(-zeta_o*w_o*Ts2)*cos(w_o*Ts2*sqrt(1-(zeta_o)^2)) == z1o + z2o;
eqn2 = exp(-2*zeta_o*w_o*Ts2) == z1o*z2o;

sol = solve([eqn1, eqn2], [z1o, z2o]);
z1o = double(sol.z1o)
z2o = double(sol.z2o)

chosen_poles_machine = ['***The chosen disc. poles (Machine) are: ', num2str(z1m(1)),' and ', num2str(z2m(1)), '***'];
chosen_poles_observer = ['***The chosen disc. poles (Observer) are: ', num2str(z1o(1)),' and ', num2str(z2o(1)), '***'];

% ------------------DIRECT DISCRETE POLES:---------------------

% % z1o = .991
% % z2o = z1o
% % 
% % p1_o = z1o + z2o; %
% % p0_o = z1o*z2o;
% 
% chosen_poles_machine = ['***The chosen disc. poles (Machine) are: ', num2str(z1m),' and ', num2str(z2m), '***'];
% chosen_poles_observer = ['***The chosen disc. poles (Observer) are: ', num2str(z1o),' and ', num2str(z2o), '***'];

%------------------------------------------------------------
%NOTE: 
% -Faster observer for modelling error rejection *
% -Slower observer for filtering high frequency noise

disp(chosen_poles_machine)
disp(chosen_poles_observer)

%% Diophantine Equ - BLANE ATTEMPT (DMC Version)
syms S2 S1 S0 r0 z
A_cl = (a2*z^2+a1*z+a0)*(z-1)*(z+r0) + (b1*z + b0)*(S2*z^2 + S1*z + S0);
A_m = (z^2-p1_m*z+p0_m);
A_o = (z^2-p1_o*z+p0_o);
A_d = A_m*A_o;
A_cl_c = fliplr(coeffs(A_cl, z)); % retreive coeficients
A_d_c = fliplr(coeffs(A_d, z));

equ1 = A_cl_c(2) == A_d_c(2); % z^3
equ2 = A_cl_c(3) == A_d_c(3); % z^2
equ3 = A_cl_c(4) == A_d_c(4); % z^1
equ4 = A_cl_c(5) == A_d_c(5); % z^0

sol = solve([equ1, equ2, equ3, equ4], [S2, S1, S0, r0]);

S2 = double(sol.S2);   
S1 = double(sol.S1);    
S0 = double(sol.S0);   
r0 = double(sol.r0);    

%% Controller TF - Gc (S/R)
%z = tf('z', Ts2);
% Gc = P + I*Ts2/(z-1) + D*(N)/(1+(N*Ts2)/(z-1)); %FOR JEREMY'S WAY

S = [S2, S1, S0];
R = [1, r0-1, -r0];

Gc = tf(S, R, Ts2)

%% Gff (T/R) and Gyr

t_o = (1-p1_m+p0_m)/(b1 + b0);
A_o = [1, -p1_o, p0_o];
A_m = [1, -p1_m, p0_m];
T = t_o*A_o;
R = [1, r0-1, -r0];

Gff = tf(T, R, Ts2)

%% MATLAB PID TUNER SUBSTITUTION

% (EULER FORWARD)
syms P I D N 
equ1 = S2 == P + (D*N);
equ2 = S1 == (P*N*Ts2) - (2*P) + (I*Ts2) - (2*D*N);
equ3 = S0 == (I*N*(Ts2^2)) - (I*Ts2) + (D*N) - (P*N*Ts2) + P;
equ4 = r0 ==(N*Ts2)-1;

sol = solve([equ4, equ1, equ2, equ3], [N, P, I, D]);

P = double(sol.P)
I = double(sol.I)
D = double(sol.D)
N = double(sol.N)

% ----- Getting b & c terms for Simulink ------------
%b = "set point weight for proportional term
%c = "set point weight for deriviative term

syms b c z
% EULEER FORWARD VERSION:
G_ff_PID = (P*b*(z-1)*(z-1+N*Ts2)) + (I*Ts2*(z-1+N*Ts2)) + (D*N*c*(z-1)^2);

% EULER BACKWARD VERSION
%G_ff_PID = (P*b*(z-1)*(z-1+N*Ts2*z)) + (I*Ts2*(z-1+N*Ts2*z)) + (D*N*c*(z-1)^2);

G_ff_c = fliplr(coeffs(G_ff_PID, z)); % retreive coeficients

equ1 = G_ff_c(1) == T(1); % z^2
equ2 = G_ff_c(2) == T(2); % z^1
equ3 = G_ff_c(3) == T(3); % z^0

sol = solve([equ1, equ3], [b, c]);

b = double(sol.b)  
c = double(sol.c)  

%% PP closed loop check (Gyr)

disp('---------- *** RESULTS - POLE PLACEMENT *** ------------')

%---------------------------ORIGINAL GYR-------------------------------

% Gyr = Gff*Gp/(1+Gc*Gp)
% poles_Gyr = pole(Gyr)
% zeros_Gyr = zero(Gyr)
% 
% Gyr = minreal(Gyr, 1e-3)
% poles_Gyr_min = pole(Gyr)
% zeros_Gyr_min = zero(Gyr)

%----------------EXPERIMENT FOR BETTER ROUNDING GYR---------------------

% Gyr = BT/(AR+BS)
syms z
A_clNew = (a2*z^2+a1*z+a0)*(z-1)*(z+r0) + (b1*z + b0)*(S2*z^2 + S1*z + S0);
A_cl_cNew = double(fliplr(coeffs(A_clNew, z))); % retreive coeficients - (AR + BS)
bT = (b2*z^2 + b1*z + b0)*(T(1)*z^2 + T(2)*z + T(3));
bT_c = double(fliplr(coeffs(bT, z)));

Gyr = tf(bT_c, A_cl_cNew, Ts2);
poles_Gyr = pole(Gyr);
zeros_Gyr = zero(Gyr);

Gyr= minreal(Gyr, 1e-3)
poles_Gyr_min = pole(Gyr)
zeros_Gyr_min = zero(Gyr)

figure(3)
pzmap(Gyr)

figure(4)
grid on
step(Gyr)

[y, t] = step(Gyr);
sserr = abs(1 - y(end));
stepinfo(Gyr)

%% AntiWindup Changes

disp('------------- *** RESULTS - PID TUNER *** --------------')

T_S2 = T(1);
T_S1 = T(2);
T_S0 = T(3);

% For I terms
c0_Gc = (S0 + S1 + S2)/(1+r0);
c0_Gff = (T_S2 + T_S1 + T_S0)/(1+r0);

% For PD Terms
d1_Gc = S2;
d1_Gff = T_S2;
d0_Gc = (r0*(S1+S2) - S0)/(1+r0);
d0_Gff = (r0*(T_S1+T_S2) - T_S0)/(1+r0); 

% PD terms:
PD_SR_num = [d1_Gc, d0_Gc];
PD_TR_num = [d1_Gff, d0_Gff];
PD_den = [1, r0];

%K_ant = 2.6; %antiwindup coefficient

%% Sensitivy Analysis

%---------------------OLD SENSITIVITY---------------------------

S_e = 1/(1 + Gp*Gc);
T_e = 1 - S_e;

%---------------------NEW SENSITIVITY (SAME RESULT)---------------------------

% S_e = (AR)/(AR+BS)
% syms z
% AR = (a2*z^2+a1*z+a0)*(z-1)*(z+r0);
% AR_c = double(fliplr(coeffs(AR, z))); % retreive coeficients - (AR + BS)
% 
% S_e = tf(AR_c, A_cl_cNew, Ts2);
% T_e = 1 - S_e;

figure(5)
bode(S_e, T_e)
legend('Sensitivity Ftn', 'Comp. Sensitivity Ftn')

figure(6)
margin(Gyr)

%% PID Tuner Parameter Method

% % Blane 2 (118.61 sec, .0152 rad/s):
% P_tun = 3.4044;
% I_tun = 0.0438;
% D_tun = 55.0446;
% N_tun = 0.0652;
% b_tun = 0.9023;
% c_tun = 0.5979;
% k_ant = 0.06;

% % Blane 4 (119.33 sec, .0202 rad/s):
% P_tun = 4.6430;
% I_tun = 0.0589;
% D_tun = 78.2077;
% N_tun = 0.0787;
% b_tun = 0.6901;
% c_tun = 0.3270;
% Kant = 0.1;

% Jeremy (119.6 sec):
P_tun  = 5.656;
I_tun  = 0.07951;
D_tun  = 6.479;
N_tun  = 0.5348;
b_tun  = 0.017;
c_tun  = 3.44E-05;
k_ant = 0.11;

% %///////////////////////////////////////////////

% Euler Fwd:
z = tf('z', Ts2);
Gc_tun = P_tun + I_tun*Ts2/(z-1) + D_tun*((N_tun)/(1+((N_tun*Ts2)/(z-1))));
Gff_tun = b_tun*P_tun + I_tun*Ts2/(z-1) + c_tun*D_tun*((N_tun)/(1+((N_tun*Ts2)/(z-1))));
Gyr_tun = Gff_tun*Gp/(1+Gc_tun*Gp);

Gyr_tun = minreal(Gyr_tun, 1e-2);

pole_tun = pole(Gyr_tun);
zero_tun = zero(Gyr_tun);

% figure(7)
% pzmap(Gyr_tun)
% 
% figure(8)
% grid on
% step(Gyr_tun)

[y, t] = step(Gyr_tun);
stepinfo(Gyr_tun);
sserr = abs(1 - y(end));


