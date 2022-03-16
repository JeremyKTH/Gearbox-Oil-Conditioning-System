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
%
T_train = getTrend(data_train);
T_test = getTrend(data_test);
T_train.InputOffset = 40;
T_train.OutputOffset = 48;
T_test.InputOffset = 40;
T_test.OutputOffset = 48;
data_train = detrend(data_train, T_train);
data_test = detrend(data_test, T_test);
%
cost_func = 'NRMSE';
%% ARMAX
opt = armaxOptions;
opt.Focus = 'prediction';
opt.SearchOptions.MaxIterations = 1000;
opt.SearchOptions.Tolerance = 1e-5;
sysARMAX = armax(data_train, [2 2 2 1], opt);
Gp = tf(sysARMAX);

[num, den] = tfdata(Gp);
b1 = num{1}(2);
b0 = num{1}(3);
a1 = den{1}(2);
a0 = den{1}(3);
B = [b1, b0];      % B
A = [1, a1, a0];  % A

% compare(data_train, sysARMAX)
%% PID Tuner - Pole Placement at 0.9904 +/- 0.0065i
P = 10.85;
I = 0.1589;
D = 168.7;
N = 2.754;
b = 0.06234;
c = 0.000279;

z = tf('z', Ts);
Gc = P + I*Ts/(z-1) + D*((N)/(1+((N*Ts)/(z-1))));
Gff = b*P + I*Ts/(z-1) + c*D*((N)/(1+((N*Ts)/(z-1))));
Gyr = Gff*Gp/(1+Gc*Gp);

Gyr = minreal(Gyr, 1e-5);

pole(Gyr);
zero(Gyr);

figure(3)
pzmap(Gyr)
figure(4)
step(Gyr)

