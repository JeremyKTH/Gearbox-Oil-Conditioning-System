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

Ts = .4; % .4 For system ID
Ts2 = .4; % After System ID
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

%% TF ESTIMATE [1 1]

opt = tfestOptions;
opt.InitializeMethod = 'all';
opt.InitialCondition = 'estimate' %Jeremy using in his code uploaded
opt.Focus = 'prediction';
opt.SearchOptions.MaxIterations = 1000; 
opt.Display = 'on';
    %opt.DistrubanceModel = 'estimate'

np = 1;          % Num of pole
nz = 1;          % Num of zero
iodelay = 0;     % **Need at 0 or cant consider for poles**

sysTF = tfest(data_train, np, nz, opt, iodelay, 'Ts', Ts);

opt = compareOptions('InitialCondition','e');

figure(1)
compare(data_test, sysTF, opt);

Gp = tf(sysTF)
Gp_cont = d2c(Gp, 'zoh'); %converting to cont.

[num, den] = tfdata(Gp); %2 poles, 1 zero
b1 = num{1}(1);     % zero z^1 (0)
b0 = num{1}(2);     % zero z^0 

a1 = den{1}(1);     % pole z^1 (1)
a0 = den{1}(2);     % pole z^0

B = [b1, b0];        % B
A = [a1, a0];  % A

%Changing sampling time for later
Gp_conv = tf(B, A, Ts2)

poles_Gp_disc = pole(Gp)
zeros_Gp_disc = zero(Gp)

figure(2)
pzmap(Gp)


%% Choose Poles (z+p_m)(z+p_o)
%--- A_m ------
a_m = .022; % continuous pole (s + a) .9(fast)---
p_m = exp(-a_m*Ts); % convert to discrete

%--- A_o ------
a_o = .018; % continuous pole (s + a)
p_o = exp(-a_o*Ts); % convert to discrete

%NOTE: 
% -Faster observer for modelling error rejection
% -Slower observer for filtering high frequency noise

%----------------------------------

chosen_disc_poles = ['**The chosen disc. poles are: ', num2str(p_m), '(machine) and ', num2str(p_o), '(observer)**'];
disp(chosen_disc_poles)


%% Diophantine Equ - BLANE ATTEMPT (DMC Version)
syms S1 S0 z
A_cl = (z+a0/a1)*(z-1) + (b1*z + b0)*(S1*z + S0);
A_d = (z - p_m)*(z - p_o);
A_cl_c = fliplr(coeffs(A_cl, z)); % retreive coeficients
A_d_c = fliplr(coeffs(A_d, z));

equ1 = A_cl_c(1) == A_d_c(1); % z^2
equ2 = A_cl_c(2) == A_d_c(2); % z^1
equ3 = A_cl_c(3) == A_d_c(3); % z^0

sol = solve([equ2, equ3], [S1, S0]);

S1 = double(sol.S1);   % 
S0 = double(sol.S0);   % 
%% Controller TF - Gc (S/R)

S = [S1, S0];
R = [1, -1];

Gc = tf(S, R, Ts)

%% Gff (T/R) and Gyr

t_o = (1 - p_m)/b0;
A_o = [1, -p_o];

T = t_o*A_o;
R = [1, -1];

Gff = tf(T, R, Ts);

%% MATLAB PID TUNER SUBSTITUTION

% (EULER FORWARD)
syms P I
equ1 = S1 == P ;
equ2 = S0 == (I*Ts - P);

sol = solve([equ1, equ2], [P, I]);

P = double(sol.P)
I = double(sol.I)

% ----- Getting b term for Simulink ------------
%b = "set point weight for proportional term

syms b z
% EULEER FORWARD VERSION:
G_ff_PI = P*b*z + (I*Ts - P*b);

% EULER BACKWARD VERSION
%G_ff_PID = (P*b*(z-1)*(z-1+N*Ts*z)) + (I*Ts*(z-1+N*Ts*z)) + (D*N*c*(z-1)^2);

G_ff_c = fliplr(coeffs(G_ff_PI, z)); % retreive coeficients

equ1 = G_ff_c(1) == T(1); % z^1
equ2 = G_ff_c(2) == T(2); % z^0

sol = solve(equ1, b);

b = double(sol)

% % Get S for Simulink
% [S, ~] = tfdata(Gc);

%% Entire closed loop check (Gyr)

disp('---------- *** RESULTS - POLE PLACEMENT *** ------------')

Gyr = Gff*Gp/(1+Gc*Gp)
% poles_Gyr = pole(Gyr)
% zeros_Gyr = zero(Gyr)

Gyr = minreal(Gyr, 1e-3)
poles_Gyr_min = pole(Gyr)
zeros_Gyr_min = zero(Gyr)

figure(3)
pzmap(Gyr)

figure(4)
grid on
step(Gyr)

[y, t] = step(Gyr);
sserr = abs(1 - y(end))
stepinfo(Gyr)

%% Sensitivy Analysis

S_e = 1/(1 + Gp*Gc);
T_e = 1 - S_e;

figure(5)
bode(S_e, T_e)
legend('Sensitivity Ftn', 'Comp. Sensitivity Ftn')

figure(6)
margin(Gyr)

%% PID Tuner Parameter Method

disp('------------- *** RESULTS - PID TUNER *** --------------')

% % Test 1 ():
% P_tun = 48.9757;
% I_tun = 1.6241;
% b_tun = .2618;
% % k_ant = 1;

% % Test 2 ():
% P_tun = 49.9724;
% I_tun = .9620;
% b_tun = .3351;
% % k_ant = 1;

% Test 3 ():
% P_tun = ;
% I_tun = ;
% b_tun = ;
% % k_ant = ;

% Jeremy ():
P_tun = 3.173;
I_tun = 0.08704;
b_tun = 0.008893;
% k_ant = ;

% %//////////////// Simulink "Tuner"  /////////////////:
% P_tun  = ;
% I_tun  = ;
% b_tun  = ;
% % k_ant = ;

% %///////////////////////////////////////////////

% Euler Forward:
z = tf('z', Ts2);
Gc_tun = P_tun  + I_tun*Ts2/(z-1);
Gff_tun = b_tun*P_tun  + I_tun *Ts2/(z-1);
Gyr_tun = Gff_tun*Gp/(1+Gc_tun*Gp);
Gyr_tun = minreal(Gyr_tun, 1e-5);

pole_tun = pole(Gyr_tun)
zero_tun = zero(Gyr_tun)

% figure(7)
% pzmap(Gyr_tun)
% 
% figure(8)
% grid on
% step(Gyr_tun)

[y, t] = step(Gyr_tun);
stepinfo(Gyr_tun)
sserr = abs(1 - y(end))

S_e = 1/(1 + Gp*Gc_tun);
T_e = 1 - S_e;

% figure(9)
% bode(S_e, T_e)
% legend('Sensitivity Ftn', 'Comp. Sensitivity Ftn')
% 
% figure(10)
% margin(Gyr)


