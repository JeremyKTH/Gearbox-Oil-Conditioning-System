%% TORQUE - ARMAX (250 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\ARMAX\ARMAX_250.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11'};

Ts = .002;
torqueRef = 250;
torqueStart = .13*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% TORQUE - ARMAX (500 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\ARMAX\ARMAX_500.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11'};

Ts = .002;
torqueRef = 500;
torqueStart = .62*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% TORQUE - ARMAX (750 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\ARMAX\ARMAX_750.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11'};

Ts = .002;
torqueRef = 750;
torqueStart = .45*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% TORQUE - ARMAX (1000 Nm)
T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\ARMAX\ARMAX_1000.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11'};

Ts = .002;
torqueRef = 1000;
torqueStart = .65*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% /////////////////// PLOT FOR TORQUE: ARMAX 250 - 1000 Nm ///////////////////////////////

% Raw data plotting variables
mixValve_Per = table2array(T_data_all(:, 1));           % OTSV1 - Mixing Valve Percent (ref.)
inGearbox_Tmp = table2array(T_data_all(:, 2));          % TV12 - Inlet Gearbox Temp
outGearbox_Tmp = table2array(T_data_all(:, 3));         % TV11 - Outlet Gearbox Temp    
inHE_Tmp = table2array(T_data_all(:, 4));               % OTGT1 - Temp into main HE (water-glycol)
outHE_Tmp = table2array(T_data_all(:, 5));              % OTGT2 - Temp out main HE (water-glycol)     
outHE_Tmp2 = table2array(T_data_all(:, 6));             % OTGT3 - Temp out main HE - further down tube (water-glycol)      
mixValveCool_Tmp = table2array(T_data_all(:, 7));       % OTGT4 - Cold Inlet Mixing Valve 
inGearbox_Flw = table2array(T_data_all(:, 8));          % QV11 - Inlet Gearbox Flow      

total_time = Ts*(length(mixValve_Per)-1);
time = (0:Ts:total_time)./60;
equil = 48;

torque = zeros(1,size(time,2));
for i=torqueStart:torqueEnd
    torque(i) = torqueRef;
end

figure(1)
clf
hold on
yyaxis left
plot(time, mixValve_Per, 'DisplayName', 'Mixing Valve [%]')
plot(time, inGearbox_Tmp, 'DisplayName', 'Oil Temp In [°C]')
plot(time, outGearbox_Tmp, 'DisplayName', 'Oil Temp Out [°C]')

%plot(time, inHE_Tmp, 'DisplayName', 'Water Temp In HE [°C]')
%plot(time, outHE_Tmp, 'DisplayName', 'Water Temp Out HE [°C]')
%plot(time, outHE_Tmp2, 'DisplayName', 'Water Temp Out HE 2 [°C]')
%plot(time, mixValveCool_Tmp, 'DisplayName', 'Cool In Mixing Valve [°C]')
% plot(time, inGearbox_Flw, 'DisplayName', 'Oil Flow [L/min.]')
yline(equil, '--m', 'DisplayName', 'Rerence [°C]')
ylabel('Parameters [misc.]')
ylim([0 75])

yyaxis right
ylabel('Torque [Nm]')
plot(time, torque, 'DisplayName', 'Torque [Nm]')
ylim([0 300])

title('Torque Disturbance')
xlabel('Time [min.]')
legend()
hold off

%% TORQUE - ARMAX (1250 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\ARMAX\ARMAX_1250.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};


%% TORQUE - SCANIA  (250 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\Scania\Scania_250.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

%% TORQUE - SCANIA  (500 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\Scania\Scania_500.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

%% TORQUE - SCANIA  (750 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\Scania\Scania_750.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

%% TORQUE - SCANIA  (1000 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\Scania\Scania_1000.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

%% TORQUE - SCANIA  (1250 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\Scania\Scania_1250.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};


%% ////////////// PLOT FOR TORQUE: ARMAX 1250 Nm, Scania 250-1250 Nm  ///////////////

% Raw data plotting variables
mixValve_Per = table2array(T_data_all(:, 1));           % OTSV1 - Mixing Valve Percent (ref.)
inGearbox_Tmp = table2array(T_data_all(:, 2));          % TV12 - Inlet Gearbox Temp
outGearbox_Tmp = table2array(T_data_all(:, 3));         % TV11 - Outlet Gearbox Temp    
inHE_Tmp = table2array(T_data_all(:, 4));               % OTGT1 - Temp into main HE (water-glycol)
outHE_Tmp = table2array(T_data_all(:, 5));              % OTGT2 - Temp out main HE (water-glycol)     
outHE_Tmp2 = table2array(T_data_all(:, 6));             % OTGT3 - Temp out main HE - further down tube (water-glycol)      
mixValveCool_Tmp = table2array(T_data_all(:, 7));       % OTGT4 - Cold Inlet Mixing Valve 

inGearbox_Flw = table2array(T_data_all(:, 8));          % QV11 - Inlet Gearbox Flow      
waterFlw = table2array(T_data_all(:, 9));                 % OTGL6 - Flow rate of oil (L/min)
coolValve_Per = table2array(T_data_all(:, 10));         % KVSV1 - Feedback for coolant valve    
inCooler_Tmp = table2array(T_data_all(:, 11));          % KVGT1 - Temp out of cooling HE  
outCooler_Tmp = table2array(T_data_all(:, 12));         % KVGT2 - Temp out of cooling HE

Ts = .002;
total_time = Ts*(length(mixValve_Per)-1);
time = (0:Ts:total_time)./60;

figure(1)
clf
hold on
plot(time, mixValve_Per, 'DisplayName', 'Mixing Valve [%]')
plot(time, inGearbox_Tmp, 'DisplayName', 'Oil Temp In [°C]')
plot(time, outGearbox_Tmp, 'DisplayName', 'Oil Temp Out [°C]')
plot(time, inHE_Tmp, 'DisplayName', 'Water Temp In HE [°C]')
plot(time, outHE_Tmp, 'DisplayName', 'Water Temp Out HE [°C]')
plot(time, outHE_Tmp2, 'DisplayName', 'Water Temp Out HE 2 [°C]')
plot(time, mixValveCool_Tmp, 'DisplayName', 'Cool In Mixing Valve [°C]')
plot(time, inGearbox_Flw, 'DisplayName', 'Oil Flow [L/min.]')
plot(time, waterFlw, 'DisplayName', 'Water Flow In HE [L/min]')
plot(time, coolValve_Per, 'DisplayName', 'Cooling Water Valve [%]')
plot(time, inCooler_Tmp, 'DisplayName', 'Cooling Water Temp In [°C]')
plot(time, outCooler_Tmp, 'DisplayName', 'Cooling Water Temp Out [°C]')
title('Torque Disturbance - 1250 Nm')
ylabel('Parameters [misc.]')
xlabel('Time (minutes)')
legend()
hold off

%% FLOW - ARMAX (350 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_350.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTGT6_DUP', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

Ts = .002;
torqueRef = 500;
torqueStart = .13*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));


%% FLOW - ARMAX (700 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_700.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTGT6_DUP', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

%% FLOW - ARMAX (1050 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_1050.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTGT6_DUP', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

%% FLOW - ARMAX (1400 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_1400.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTGT6_DUP', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

%% FLOW - ARMAX (2100 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_2100.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'KVSV1', 'KVGT1', 'KVGT2'};

%% FLOW - ARMAX (2450 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_2450.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTGT6_DUP', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

%% FLOW - ARMAX (2800 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_2800.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11', 'KVSV1', 'KVGT1', 'KVGT2'};

%% FLOW - ARMAX (3150 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_3150.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTGT6_DUP', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

%% FLOW - ARMAX (3500 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_3500.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'QV11', 'OTGT6_DUP', 'KVSV1', 'KVGT1', 'KVGT2', 'OTGP1'};

%% ////////////// PLOT FOR FLOW: ARMAX ///////////////

% Raw data plotting variables
mixValve_Per = table2array(T_data_all(:, 1));           % OTSV1 - Mixing Valve Percent (ref.)
inGearbox_Tmp = table2array(T_data_all(:, 2));          % TV12 - Inlet Gearbox Temp
outGearbox_Tmp = table2array(T_data_all(:, 3));         % TV11 - Outlet Gearbox Temp    
inHE_Tmp = table2array(T_data_all(:, 4));               % OTGT1 - Temp into main HE (water-glycol)
outHE_Tmp = table2array(T_data_all(:, 5));              % OTGT2 - Temp out main HE (water-glycol)     
outHE_Tmp2 = table2array(T_data_all(:, 6));             % OTGT3 - Temp out main HE - further down tube (water-glycol)      
mixValveCool_Tmp = table2array(T_data_all(:, 7));       % OTGT4 - Cold Inlet Mixing Valve 

mixValveHot_Tmp = table2array(T_data_all(:, 8));        % OTGT6 - Hot Inlet Mixing Valve 
inGearbox_Flw = table2array(T_data_all(:, 9));          % QV11 - Inlet Gearbox Flow
mixValveHot_Tmp2 = table2array(T_data_all(:, 10));       % OTGT6 - Hot Inlet Mixing Valve (DUPLICATE??)
% waterFlw = table2array(T_data_all(:, 10));            % OTGL6 - Flow rate of oil (L/min)
coolValve_Per = table2array(T_data_all(:, 11));         % KVSV1 - Feedback for coolant valve    
inCooler_Tmp = table2array(T_data_all(:, 12));          % KVGT1 - Temp out of cooling HE  
outCooler_Tmp = table2array(T_data_all(:, 13));         % KVGT2 - Temp out of cooling HE
inHE_Pre = table2array(T_data_all(:, 14));              % OTGP1 - Pressure out of mixing valve


total_time = Ts*(length(mixValve_Per)-1);
time = (0:Ts:total_time)./60;
equil = 48;

torque = zeros(1,size(time,2));
for i=torqueStart:torqueEnd
    torque(i) = torqueRef;
end

figure(1)
clf
hold on
plot(time, mixValve_Per, 'DisplayName', 'Mixing Valve [%]') % OTSV1 - Mixing Valve Percent (ref.)
plot(time, inGearbox_Tmp, 'DisplayName', 'Oil Temp In [°C]') % TV12 - Inlet Gearbox Temp
plot(time, outGearbox_Tmp, 'DisplayName', 'Oil Temp Out [°C]') % TV11 - Outlet Gearbox Temp 
plot(time, inHE_Tmp, 'DisplayName', 'Water Temp In HE [°C]')  % OTGT1 - Temp into main HE (water-glycol)
plot(time, outHE_Tmp, 'DisplayName', 'Water Temp Out HE [°C]') % OTGT2 - Temp out main HE (water-glycol) 
plot(time, outHE_Tmp2, 'DisplayName', 'Water Temp Out HE 2 [°C]') % OTGT3 - Temp out main HE - further down tube (water-glycol)  
plot(time, mixValveCool_Tmp, 'DisplayName', 'Cool In Mixing Valve [°C]') % OTGT4 - Cold Inlet Mixing Valve 

plot(time, mixValveHot_Tmp, 'DisplayName', 'Hot In Mixing Valve [°C]') % OTGT6 - Hot Inlet Mixing Valve 
plot(time, inGearbox_Flw, 'DisplayName', 'Oil Flow [L/min.]')
plot(time, mixValveHot_Tmp2, 'DisplayName', '**Hot In Mixing Valve** [°C]')
% plot(time, waterFlw, 'DisplayName', 'Water Flow In HE [L/min]') % OTGL6 - Flow rate of oil (L/min)
plot(time, coolValve_Per, 'DisplayName', 'Cooling Water Valve [%]')
plot(time, inCooler_Tmp, 'DisplayName', 'Cooling Water Temp In [°C]')
plot(time, outCooler_Tmp, 'DisplayName', 'Cooling Water Temp Out [°C]')
plot(time, inHE_Pre, 'DisplayName', 'Cooling Water Temp Out [Pa]')

plot(time, torque, 'DisplayName', 'Torque [Nm]')
yline(equil, '--m', 'DisplayName', 'Reference [L/min.]')
ylim([0 100])
title('Flow Disturbance - 350 RPM')
xlabel('Time (minutes)')
legend()
hold off

%% FLOW - SCANIA

% --- 350 RPM ---

% --- 700 RPM ---

% --- 1050 RPM ---

% --- 1400 RPM ---

% --- 2100 RPM ---

% --- 2450 RPM ---

% --- 2800 RPM ---

% --- 3150 RPM ---

% --- 3500 RPM ---

