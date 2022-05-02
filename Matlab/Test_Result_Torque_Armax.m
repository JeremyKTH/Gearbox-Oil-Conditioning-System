%% TORQUE - ARMAX (250 Nm)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Torque\ARMAX\ARMAX_250.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_250  = mean(T_data_all.TV12)
std_ARMAX_250   = std(T_data_all.TV12)
nrmse_ARMAX_250 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_250

%% TORQUE - ARMAX (500 Nm)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Torque\ARMAX\ARMAX_500.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_500  = mean(T_data_all.TV12)
std_ARMAX_500   = std(T_data_all.TV12)
nrmse_ARMAX_500 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_500

%% TORQUE - ARMAX (750 Nm)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Torque\ARMAX\ARMAX_750.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_750  = mean(T_data_all.TV12)
std_ARMAX_750   = std(T_data_all.TV12)
nrmse_ARMAX_750 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_750

%% TORQUE - ARMAX (1000 Nm)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Torque\ARMAX\ARMAX_1000.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_1000  = mean(T_data_all.TV12)
std_ARMAX_1000   = std(T_data_all.TV12)
nrmse_ARMAX_1000 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_1000

%% TORQUE - ARMAX (1250 Nm)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Torque\ARMAX\ARMAX_1250.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11', 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_1250  = mean(T_data_all.TV12)
std_ARMAX_1250   = std(T_data_all.TV12)
nrmse_ARMAX_1250 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_1250

%% TORQUE - TOTAL
mean_ARMAX_total = (mean_ARMAX_250+mean_ARMAX_500+mean_ARMAX_750+mean_ARMAX_1000+mean_ARMAX_1250)/5
std_ARMAX_total = (std_ARMAX_250+std_ARMAX_500+std_ARMAX_750+std_ARMAX_1000+std_ARMAX_1250)/5
nrmse_ARMAX_total = (nrmse_ARMAX_250+nrmse_ARMAX_500+nrmse_ARMAX_750+nrmse_ARMAX_1000+nrmse_ARMAX_1250)/5