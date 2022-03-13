T_train = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\Test_Data\Training_Data\D_48404_train.csv');
T_train.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1'};
T_test = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\Test_Data\Testing_Data\D_48404_test.csv');
T_test.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1'};
u_train = table2array(T_train(:, 1));   % OTSV1
y_train = table2array(T_train(:, 2));   % TV12
u_test =  table2array(T_test(:, 1));   % OTSV1
y_test =  table2array(T_test(:, 2));   % TV12
Ts = 0.4;  % second
% Training Data
data_train = iddata(y_train, u_train, Ts);
data_train.Name = 'Data_Train';
data_train.TimeUnit = 'seconds';
data_train.InputName = 'OTSV1';   data_train.InputUnit = 'Percentage';
data_train.OutputName = 'TV12';   data_train.OutputUnit = 'Celsius';
% Testing Data
data_test = iddata(y_test, u_test, Ts);
data_test.Name = 'Data_Test';
data_test.TimeUnit = 'seconds';
data_test.InputName = 'OTSV1';   data_test.InputUnit = 'Percentage';
data_test.OutputName = 'TV12';   data_test.OutputUnit = 'Celsius';
data_train = detrend(data_train);
data_test = detrend(data_test);
%
cost_func = 'NRMSE';
%% TF
opt = tfestOptions;
opt.InitializeMethod = 'all';
opt.SearchOptions.MaxIterations = 1000;
iodelay = 0.0;   % In/Out delay
Gp = tfest(data_train, 1, 1, opt, iodelay, 'Ts', Ts);

[num, den] = tfdata(Gp);
b = num{1}(2);
a = den{1}(2);
B = [b];     % B
A = [1, a];  % A
%% Best +/-40% TF Model
% b = 0.001273;
% a = -0.9965;
% B = [b];         % B
% A = [1, a];  % A
% Ts = 0.4;
% 
% Gp = tf(B, A, Ts);
% 
% poles = pole(Gp);  % 0.9965
% zeros = zero(Gp);  % No zeros
%pzmap(Gp)
%step(Gp)
%% Choose Poles (w_m) (w_o)
%--- A_m ------
w_m = 1;
%--- A_o ------
w_o = 1;

p1 = exp(-w_m*Ts);
p0 = exp(-w_o*Ts);

fprintf('The original disc. poles are: 0.9965 \n');
fprintf('The chosen cont. poles are: %d and %d \n', -p1, -p0);
chosen_disc_poles = ['The chosen disc. poles are: ', num2str(p1), ' and ', num2str(p0)];
disp(chosen_disc_poles)

%% Diophantine Eqn A_cl = AR + BS
syms P I z
A_cl = (z+a)*((z-1)*P + I*Ts);
A_d = (z-p1)*(z-p0);
A_cl_c = fliplr(coeffs(A_cl, z)); % retreive coeficients
A_d_c = fliplr(coeffs(A_d, z));   % retreive coeficients

equ1 = A_cl_c(2) == A_d_c(2); % z^1
equ2 = A_cl_c(3) == A_d_c(3); % z^0

sol = solve([equ1, equ2], [P, I]);

P = double(sol.P)   % 0.8929
I = double(sol.I)   % 1.1049

%% Gc
z = tf('z', Ts);
Gc = P + I*Ts/(z-1);

%% Gff and Gyr
t_o = (1-p1)/(b);
A_o = [1, -p0];
T = t_o*A_o;
R = [1, -1];
Gff = tf(T, R, Ts);

%% Entire closed loop check
Gyr = Gff*Gp/(1+Gc*Gp)
pole(Gyr)
zero(Gyr)

figure(1)
pzmap(Gyr)

figure(2)
step(Gyr)
%% The following parts are just for Simulink Usage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve b & c for Simulink
% syms b z
% G_ff =  (z-1)*b*P+I*Ts;
% G_ff_c = fliplr(coeffs(G_ff, z)); % retreive coeficients
% 
% equ1 = G_ff_c(1) == T(1); % z^1
% equ2 = G_ff_c(2) == T(2); % z^0
% 
% 
% sol = solve([equ1], [b]);
% 
% b = double(sol.b)   % 0.4684
%b = T(1)/P
b = -(T(2)-I*Ts)/P
%% Get S for Simulink
[S, ~] = tfdata(Gc);

%% PID Tuner
P = 4.9013;
I = 0.042267;
b = 1;