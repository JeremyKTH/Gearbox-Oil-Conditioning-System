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
%% ARX
opt = arxOptions;
opt.Focus = 'prediction';
sysARX = arx(data_train, [3, 1, 1])  %------------------------------------
y_arx_ref = data_test.y;
y_arx_est = sim(sysARX, data_test(:, [], :));
fit = goodnessOfFit(y_arx_est.y, y_arx_ref, cost_func)
FPE = sysARX.Report.Fit.FPE
figure(1)
resid(data_test,sysARX)
figure(2)
compare(data_test, sysARX)
%% ARMAX
opt = armaxOptions;
opt.Focus = 'prediction';
opt.SearchOptions.MaxIterations = 1000;
opt.SearchOptions.Tolerance = 1e-5;
sysARMAX = armax(data_train, [2 1 1 1], opt) %----------------------------
y_armax_ref = data_test.y;
y_armax_est = sim(sysARMAX, data_test(:, [], :));
fit = goodnessOfFit(y_armax_est.y, y_armax_ref, cost_func)
FPE = sysARMAX.Report.Fit.FPE
figure(1)
resid(data_test,sysARMAX)
figure(2)
compare(data_test, sysARMAX)
%% BJ
opt = bjOptions;
opt.Focus = 'prediction';
opt.SearchOptions.MaxIterations = 1000;
opt.SearchOptions.Tolerance = 1e-5;
sysBJ = bj(data_train, [1 1 2 1 1], opt, 'IODelay', 20); % ----------------
y_bj_est = sim(sysBJ, data_test(:, [], :));
y_bj_ref = data_test.y;
fit = goodnessOfFit(y_bj_est.y, y_bj_ref, cost_func)
FPE = sysBJ.Report.Fit.FPE
figure(1)
resid(data_test,sysBJ)
figure(2)
compare(data_test, sysBJ)
%% TF
opt = tfestOptions;
opt.InitializeMethod = 'all';
opt.SearchOptions.MaxIterations = 1000;
iodelay = 20.0;   % In/Out delay
sysTF = tfest(data_train, 1, 0, opt, iodelay, 'Ts', Ts) % ----------------
y_tf_ref = data_test.y;
y_tf_est = sim(sysTF, data_test(:, [], :));
fit = goodnessOfFit(y_tf_est.y, y_tf_ref, cost_func)
FPE = sysTF.Report.Fit.FPE
figure(1)
resid(data_test,sysTF)
figure(2)
compare(data_test,sysTF)
%% SS
CoarseOrderSearch = [1, 2, 3, 4];
sysSS = n4sid(data_train, CoarseOrderSearch)
%%
y_ss_est = sim(sysSS, data_test(:, [], :)); % -----------------------------
y_ss_ref = data_test.y;
fit = goodnessOfFit(y_ss_est.y, y_ss_ref, cost_func)
FPE = sysSS.Report.Fit.FPE
figure(1)
resid(data_test,sysSS)
figure(2)
compare(data_test,sysSS)