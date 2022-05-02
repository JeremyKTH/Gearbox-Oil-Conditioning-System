%% Flow - ARMAX (350 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\ARMAX\ARMAX_350.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_350  = mean(T_data_all.TV12)
std_ARMAX_350   = std(T_data_all.TV12)
nrmse_ARMAX_350 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_350

%% Flow - ARMAX (700 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\ARMAX\ARMAX_700.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_700  = mean(T_data_all.TV12)
std_ARMAX_700   = std(T_data_all.TV12)
nrmse_ARMAX_700 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_700

%% Flow - ARMAX (1050 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\ARMAX\ARMAX_1050.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_1050  = mean(T_data_all.TV12)
std_ARMAX_1050   = std(T_data_all.TV12)
nrmse_ARMAX_1050 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_1050

%% Flow - ARMAX (1400 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\ARMAX\ARMAX_1400.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_1400  = mean(T_data_all.TV12)
std_ARMAX_1400   = std(T_data_all.TV12)
nrmse_ARMAX_1400 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_1400

%% Flow - ARMAX (2100 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\ARMAX\ARMAX_2100.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11', 'KVSV1', 'KVGT1', 'KVGT2'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_2100  = mean(T_data_all.TV12)
std_ARMAX_2100   = std(T_data_all.TV12)
nrmse_ARMAX_2100 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_2100

%% Flow - ARMAX (2450 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\ARMAX\ARMAX_2450.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_2450  = mean(T_data_all.TV12)
std_ARMAX_2450   = std(T_data_all.TV12)
nrmse_ARMAX_2450 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_2450

%% Flow - ARMAX (2800 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\ARMAX\ARMAX_2800.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_2800  = mean(T_data_all.TV12)
std_ARMAX_2800   = std(T_data_all.TV12)
nrmse_ARMAX_2800 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_2800

%% Flow - ARMAX (3150 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\ARMAX\ARMAX_3150.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_3150  = mean(T_data_all.TV12)
std_ARMAX_3150   = std(T_data_all.TV12)
nrmse_ARMAX_3150 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_3150

%% Flow - ARMAX (3500 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\ARMAX\ARMAX_3500.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_ARMAX_3500  = mean(T_data_all.TV12)
std_ARMAX_3500   = std(T_data_all.TV12)
nrmse_ARMAX_3500 = sqrt(mse(T_data_all.TV12, ref))/mean_ARMAX_3500

%% Flow - TOTAL
mean_ARMAX_total = (mean_ARMAX_350+mean_ARMAX_700+mean_ARMAX_1050+mean_ARMAX_1400+mean_ARMAX_2100+mean_ARMAX_2450+mean_ARMAX_2800+mean_ARMAX_3150+mean_ARMAX_3500)/9
std_ARMAX_total = (std_ARMAX_350+std_ARMAX_700+std_ARMAX_1050+std_ARMAX_1400+std_ARMAX_2100+std_ARMAX_2450+std_ARMAX_2800+std_ARMAX_3150+std_ARMAX_3500)/9
nrmse_ARMAX_total = (nrmse_ARMAX_350+nrmse_ARMAX_700+nrmse_ARMAX_1050+nrmse_ARMAX_1400+nrmse_ARMAX_2100+nrmse_ARMAX_2450+nrmse_ARMAX_2800+nrmse_ARMAX_3150+nrmse_ARMAX_3500)/9