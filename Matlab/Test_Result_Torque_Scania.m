%% TORQUE - Scania (250 Nm)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Torque\Scania\Scania_250.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11', 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_250  = mean(T_data_all.TV12)
std_Scania_250   = std(T_data_all.TV12)
nrmse_Scania_250 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_250

%% TORQUE - Scania (500 Nm)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Torque\Scania\Scania_500.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11', 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_500  = mean(T_data_all.TV12)
std_Scania_500   = std(T_data_all.TV12)
nrmse_Scania_500 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_500

%% TORQUE - Scania (750 Nm)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Torque\Scania\Scania_750.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11', 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_750  = mean(T_data_all.TV12)
std_Scania_750   = std(T_data_all.TV12)
nrmse_Scania_750 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_750

%% TORQUE - Scania (1000 Nm)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Torque\Scania\Scania_1000.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11', 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_1000  = mean(T_data_all.TV12)
std_Scania_1000   = std(T_data_all.TV12)
nrmse_Scania_1000 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_1000

%% TORQUE - Scania (1250 Nm)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Torque\Scania\Scania_1250.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11', 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_1250  = mean(T_data_all.TV12)
std_Scania_1250   = std(T_data_all.TV12)
nrmse_Scania_1250 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_1250

%% TORQUE - TOTAL
mean_Scania_total = (mean_Scania_250+mean_Scania_500+mean_Scania_750+mean_Scania_1000+mean_Scania_1250)/5
std_Scania_total = (std_Scania_250+std_Scania_500+std_Scania_750+std_Scania_1000+std_Scania_1250)/5
nrmse_Scania_total = (nrmse_Scania_250+nrmse_Scania_500+nrmse_Scania_750+nrmse_Scania_1000+nrmse_Scania_1250)/5