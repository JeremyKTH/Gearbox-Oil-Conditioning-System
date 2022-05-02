%% Flow - Scania (350 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\Scania\Scania_350.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_350  = mean(T_data_all.TV12)
std_Scania_350   = std(T_data_all.TV12)
nrmse_Scania_350 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_350

%% Flow - Scania (700 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\Scania\Scania_700.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_700  = mean(T_data_all.TV12)
std_Scania_700   = std(T_data_all.TV12)
nrmse_Scania_700 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_700

%% Flow - Scania (1050 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\Scania\Scania_1050.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_1050  = mean(T_data_all.TV12)
std_Scania_1050   = std(T_data_all.TV12)
nrmse_Scania_1050 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_1050

%% Flow - Scania (1400 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\Scania\Scania_1400.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_1400  = mean(T_data_all.TV12)
std_Scania_1400   = std(T_data_all.TV12)
nrmse_Scania_1400 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_1400

%% Flow - Scania (2100 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\Scania\Scania_2100.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_2100  = mean(T_data_all.TV12)
std_Scania_2100   = std(T_data_all.TV12)
nrmse_Scania_2100 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_2100

%% Flow - Scania (2450 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\Scania\Scania_2450.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_2450  = mean(T_data_all.TV12)
std_Scania_2450   = std(T_data_all.TV12)
nrmse_Scania_2450 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_2450

%% Flow - Scania (2800 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\Scania\Scania_2800.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_2800  = mean(T_data_all.TV12)
std_Scania_2800   = std(T_data_all.TV12)
nrmse_Scania_2800 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_2800

%% Flow - Scania (3150 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\Scania\Scania_3150.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_3150  = mean(T_data_all.TV12)
std_Scania_3150   = std(T_data_all.TV12)
nrmse_Scania_3150 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_3150

%% Flow - Scania (3500 L/min)
T_data_all = readtable('C:\Users\jerem\Documents\Python Scripts\Scania\RQ2_Data\Clean_Data\Flow\Scania\Scania_3500.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTSV1_FB_LREAL', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

ref = zeros(length(T_data_all.TV12), 1);
ref(1:end)= 48;

mean_Scania_3500  = mean(T_data_all.TV12)
std_Scania_3500   = std(T_data_all.TV12)
nrmse_Scania_3500 = sqrt(mse(T_data_all.TV12, ref))/mean_Scania_3500

%% Flow - TOTAL
mean_Scania_total = (mean_Scania_350+mean_Scania_700+mean_Scania_1050+mean_Scania_1400+mean_Scania_2100+mean_Scania_2450+mean_Scania_2800+mean_Scania_3150+mean_Scania_3500)/9
std_Scania_total = (std_Scania_350+std_Scania_700+std_Scania_1050+std_Scania_1400+std_Scania_2100+std_Scania_2450+std_Scania_2800+std_Scania_3150+std_Scania_3500)/9
nrmse_Scania_total = (nrmse_Scania_350+nrmse_Scania_700+nrmse_Scania_1050+nrmse_Scania_1400+nrmse_Scania_2100+nrmse_Scania_2450+nrmse_Scania_2800+nrmse_Scania_3150+nrmse_Scania_3500)/9