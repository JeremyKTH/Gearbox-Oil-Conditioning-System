%% TORQUE - ARMAX (250 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\ARMAX\ARMAX_250.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11'};

Ts = .002;
ctrl = '(ARMAX Derived)';
torqueRef = 250;
torqueStart = .13*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% TORQUE - ARMAX (500 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\ARMAX\ARMAX_500.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11'};

Ts = .002;
ctrl = '(ARMAX Derived)';
torqueRef = 500;
torqueStart = .62*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% TORQUE - ARMAX (750 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\ARMAX\ARMAX_750.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11'};

Ts = .002;
ctrl = '(ARMAX Derived)';
torqueRef = 750;
torqueStart = .45*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% TORQUE - ARMAX (1000 Nm)
T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\ARMAX\ARMAX_1000.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11'};

Ts = .002;
ctrl = '(ARMAX Derived)';
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
plot(time, mixValve_Per, '-', 'Color',[0 0.4470 0.7410], 'DisplayName', 'Mixing Valve [%]', 'LineWidth', 2)
plot(time, inGearbox_Tmp, '-', 'Color',[0.8500 0.3250 0.0980], 'DisplayName', 'Oil Temp In [°C]', 'LineWidth', 2)
plot(time, outGearbox_Tmp, '-', 'Color',[0.9290 0.6940 0.1250], 'DisplayName', 'Oil Temp Out [°C]', 'LineWidth', 2)
% plot(time, inHE_Tmp, 'DisplayName', 'Water Temp In HE [°C]')
% plot(time, outHE_Tmp, 'DisplayName', 'Water Temp Out HE [°C]')
% plot(time, outHE_Tmp2, 'DisplayName', 'Water Temp Out HE 2 [°C]')
% plot(time, mixValveCool_Tmp, 'DisplayName', 'Cool In Mixing Valve [°C]')
% plot(time, inGearbox_Flw, 'DisplayName', 'Oil Flow [L/min.]')
yline(equil, '--k', 'DisplayName', 'Reference - 48 °C', 'LineWidth', 2)

ylabel('Temperature [°C] and Percentage [%]')
ylim([(round(min(mixValve_Per),-1)-5) (round(max(outGearbox_Tmp), -1)+5)]);  
yticks((round(min(mixValve_Per),-1)-5):5:(round(max(outGearbox_Tmp), -1)+5));

yyaxis right
plot(time, torque, ':', 'Color',[0.4940 0.1840 0.5560], 'DisplayName', 'Torque [Nm]', 'LineWidth', 2)
ylabel('Torque [Nm]')
ylim([0 torqueRef+25]);


title(['Torque Disturbance - ', num2str(torqueRef), ' Nm ', ctrl])
xlabel('Time [min.]')
xlim([0 max(time)])
ax = gca;
ax.FontSize = 18;
ax.YColor = [0.4940 0.1840 0.5560];
legend()
hold off

%% TORQUE - ARMAX (1250 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\ARMAX\ARMAX_1250.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

Ts = .002;
ctrl = '(ARMAX Derived)';
torqueRef = 1250;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = .65*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% TORQUE - SCANIA  (250 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\Scania\Scania_250.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

Ts = .002;
ctrl = '(Manually Tuned)';
torqueRef = 250;
simStart = double(5*(60*(1/Ts)));
torqueStart = .16*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% TORQUE - SCANIA  (500 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\Scania\Scania_500.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

Ts = .002;
ctrl = '(Manually Tuned)';
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = .43*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% TORQUE - SCANIA  (750 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\Scania\Scania_750.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

Ts = .002;
ctrl = '(Manually Tuned)';
torqueRef = 750;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = .35*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% TORQUE - SCANIA  (1000 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\Scania\Scania_1000.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

Ts = .002;
ctrl = '(Manually Tuned)';
torqueRef = 1000;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = .35*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% TORQUE - SCANIA  (1250 Nm)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Torque\Scania\Scania_1250.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'QV11' , 'OTGL6', 'KVSV1', 'KVGT1', 'KVGT2'};

Ts = .002;
ctrl = '(Manually Tuned)';
torqueRef = 1250;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = .50*(60*(1/Ts)); % .13 min start
torqueEnd = torqueStart + (5*60*(1/Ts));

%% ////////////// PLOT FOR TORQUE: ARMAX 1250 Nm, Scania 250-1250 Nm  ///////////////

% Raw data plotting variables
mixValve_Per = table2array(T_data_all(simStart:end, 1));           % OTSV1 - Mixing Valve Percent (ref.)
inGearbox_Tmp = table2array(T_data_all(simStart:end, 2));          % TV12 - Inlet Gearbox Temp
outGearbox_Tmp = table2array(T_data_all(simStart:end, 3));         % TV11 - Outlet Gearbox Temp    
inHE_Tmp = table2array(T_data_all(simStart:end, 4));               % OTGT1 - Temp into main HE (water-glycol)
outHE_Tmp = table2array(T_data_all(simStart:end, 5));              % OTGT2 - Temp out main HE (water-glycol)     
outHE_Tmp2 = table2array(T_data_all(simStart:end, 6));             % OTGT3 - Temp out main HE - further down tube (water-glycol)      
mixValveCool_Tmp = table2array(T_data_all(simStart:end, 7));       % OTGT4 - Cold Inlet Mixing Valve 

inGearbox_Flw = table2array(T_data_all(simStart:end, 8));          % QV11 - Inlet Gearbox Flow      
waterFlw = table2array(T_data_all(simStart:end, 9));                 % OTGL6 - Flow rate of oil (L/min)
coolValve_Per = table2array(T_data_all(simStart:end, 10));         % KVSV1 - Feedback for coolant valve    
inCooler_Tmp = table2array(T_data_all(simStart:end, 11));          % KVGT1 - Temp out of cooling HE  
outCooler_Tmp = table2array(T_data_all(simStart:end, 12));         % KVGT2 - Temp out of cooling HE

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
plot(time, mixValve_Per, '-', 'Color',[0 0.4470 0.7410], 'DisplayName', 'Mixing Valve [%]', 'LineWidth', 2)
plot(time, inGearbox_Tmp, '-', 'Color',[0.8500 0.3250 0.0980], 'DisplayName', 'Oil Temp In [°C]', 'LineWidth', 2)
plot(time, outGearbox_Tmp, '-', 'Color',[0.9290 0.6940 0.1250], 'DisplayName', 'Oil Temp Out [°C]', 'LineWidth', 2)
% plot(time, inHE_Tmp, 'DisplayName', 'Water Temp In HE [°C]')
% plot(time, outHE_Tmp, 'DisplayName', 'Water Temp Out HE [°C]')
% plot(time, outHE_Tmp2, 'DisplayName', 'Water Temp Out HE 2 [°C]')
% plot(time, mixValveCool_Tmp, 'DisplayName', 'Cool In Mixing Valve [°C]')
% plot(time, inGearbox_Flw, 'DisplayName', 'Oil Flow [L/min.]')
% plot(time, waterFlw, 'DisplayName', 'Water Flow In HE [L/min]')
% plot(time, coolValve_Per, 'DisplayName', 'Cooling Water Valve [%]')
% plot(time, inCooler_Tmp, 'DisplayName', 'Cooling Water Temp In [°C]')
% plot(time, outCooler_Tmp, 'DisplayName', 'Cooling Water Temp Out [°C]')
yline(equil, '--k', 'DisplayName', 'Reference - 48 °C', 'LineWidth', 2)

ylabel('Temperature [°C] and Percentage [%]')
ylim([(round(min(mixValve_Per),-1)-5) (round(max(outGearbox_Tmp), -1)+5)]);  
yticks((round(min(mixValve_Per),-1)-5):5:(round(max(outGearbox_Tmp), -1)+5));

yyaxis right
plot(time, torque, ':', 'Color',[0.4940 0.1840 0.5560], 'DisplayName', 'Torque [Nm]', 'LineWidth', 2)
ylabel('Torque [Nm]')
ylim([0 torqueRef+25]);


title(['Torque Disturbance - ', num2str(torqueRef), ' Nm ', ctrl])
xlabel('Time [min.]')
xlim([0 max(time)])
ax = gca;
ax.FontSize = 18;
ax.YColor = [0.4940 0.1840 0.5560];
legend()
hold off

%% FLOW - ARMAX (350 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_350.csv');
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(ARMAX Derived)';
flowRef = 350;
torqueRef = 500;
torqueStart = double(.556*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));


%% FLOW - ARMAX (700 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_700.csv');
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(ARMAX Derived)';
flowRef = 700;
torqueRef = 500;
simStart = double(4.5*(60*(1/Ts)));
torqueStart = double(.68*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% FLOW - ARMAX (1050 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_1050.csv');
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(ARMAX Derived)';
flowRef = 1050;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(1.56*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% FLOW - ARMAX (1400 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_1400.csv');
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(ARMAX Derived)';
flowRef = 1400;
torqueRef = 500;
simStart = double(7.2*(60*(1/Ts)));
torqueStart = double(1.73*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% FLOW - ARMAX (2100 RPM) ** SPECIAL CASE IMPORT **

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_2100.csv');
T_data_all.Properties.VariableNames = {'X1', 'TV12', 'X2', 'TV11', 'X3' 'QV11', 'X4', 'OTGT4', 'X5', 'OTGT2', 'X6', 'OTGT1', 'X7', 'OTGT3', 'X8', 'OTSV1', 'X9', 'OTGL6', 'X10', 'KVSV1', 'X11', 'KVGT1', 'X12', 'KVGT2'};

Ts = .002;
ctrl = '(ARMAX Derived)';
flowRef = 2100;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(.5*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));


%% FLOW - ARMAX (2450 RPM) 

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_2450.csv');
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(ARMAX Derived)';
flowRef = 2450;
torqueRef = 500;
simStart = double(.1*(60*(1/Ts)));
torqueStart = double(.35*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% FLOW - ARMAX (2800 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_2800.csv');
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(ARMAX Derived)';
flowRef = 2800;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(.4*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% FLOW - ARMAX (3150 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_3150.csv');
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(ARMAX Derived)';
flowRef = 3150;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(.59*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% FLOW - ARMAX (3500 RPM)

T_data_all = readtable('H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\ARMAX\ARMAX_3500.csv');
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(ARMAX Derived)';
flowRef = 3500;
torqueRef = 500;
simStart = double(.5*(60*(1/Ts)));
torqueStart = double(.86*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% ////////////// PLOT FOR FLOW: ARMAX ///////////////

% Raw data plotting variables (For all except 2100)
inCooler_Tmp = table2array(T_data_all(simStart:end, 2));           % KVGT1 - Temp out of cooling HE  
outCooler_Tmp = table2array(T_data_all(simStart:end, 4));          % KVGT2 - Temp out of cooling HE
inHE_Tmp = table2array(T_data_all(simStart:end, 6));               % OTGT1 - Temp into main HE (water-glycol)
outHE_Tmp = table2array(T_data_all(simStart:end, 8));              % OTGT2 - Temp out main HE (water-glycol) 
outHE_Tmp2 = table2array(T_data_all(simStart:end, 10));             % OTGT3 - Temp out main HE - further down tube (water-glycol)      
mixValveCool_Tmp = table2array(T_data_all(simStart:end, 12));       % OTGT4 - Cold Inlet Mixing Valve 
mixValveHot_Tmp = table2array(T_data_all(simStart:end, 14));        % OTGT6 - Hot Inlet Mixing Valve 
mixValve_Per = table2array(T_data_all(simStart:end, 16));           % OTSV1 - Mixing Valve Percent (ref.)
inGearbox_Flw = table2array(T_data_all(simStart:end, 18));          % QV11 - Inlet Gearbox Flow
outGearbox_Tmp = table2array(T_data_all(simStart:end, 20));        % TV11 - Outlet Gearbox Temp  
inGearbox_Tmp = table2array(T_data_all(simStart:end, 22));         % TV12 - Inlet Gearbox Temp
ref_Tmp = table2array(T_data_all(simStart:end, 24));               % TV_Cmd - Reference TV12 Temp
mixValve_REAL = table2array(T_data_all(simStart:end, 26));         % OTSV1_FB_LREAL - REAL Mixing Valve Percent
inHE_Pre = table2array(T_data_all(simStart:end, 28));              % OTGP1 - Pressure out of mixing valve
coolValve_Per = table2array(T_data_all(simStart:end, 30));         % KVSV1 - Coolant Valve Percent    

% % Raw data plotting variables *** 2100 RPM SPECIAL CASE ONLY ***
% inGearbox_Tmp = table2array(T_data_all(simStart:end, 2));          % TV12 - Inlet Gearbox Temp
% outGearbox_Tmp = table2array(T_data_all(simStart:end, 4));         % TV11 - Outlet Gearbox Temp
% inGearbox_Flw = table2array(T_data_all(simStart:end, 6));          % QV11 - Inlet Gearbox Flow
% mixValveCool_Tmp = table2array(T_data_all(simStart:end, 8));       % OTGT4 - Cold Inlet Mixing Valve 
% outHE_Tmp = table2array(T_data_all(simStart:end, 10));               % OTGT2 - Temp out main HE (water-glycol) 
% inHE_Tmp = table2array(T_data_all(simStart:end, 12));                % OTGT1 - Temp into main HE (water-glycol)
% outHE_Tmp2 = table2array(T_data_all(simStart:end, 14));             % OTGT3 - Temp out main HE - further down tube (water-glycol) 
% mixValve_Per = table2array(T_data_all(simStart:end, 16));           % OTSV1 - Mixing Valve Percent (ref.)
% waterFlw = table2array(T_data_all(:, 18));                           % OTGL6 - Flow rate of oil (L/min)
% coolValve_Per = table2array(T_data_all(simStart:end, 20));          % KVSV1 - Coolant Valve Percent    
% inCooler_Tmp = table2array(T_data_all(simStart:end, 22));            % KVGT1 - Temp out of cooling HE  
% outCooler_Tmp = table2array(T_data_all(simStart:end, 24));           % KVGT2 - Temp out of cooling HE

total_time = Ts*(length(mixValve_Per)-1);
time = (0:Ts:total_time)./60;
equil = 48;

% torque = zeros(1,size(time,2));
% for i=torqueStart:torqueEnd
%     torque(i) = torqueRef;
% end

flow = zeros(1,size(time,2))+1750;
for i=torqueStart:torqueEnd
    flow(i) = flowRef;
end

figure(1)
clf
hold on
ax = gca;
ax.FontSize = 18;
yyaxis left
% plot(time, inCooler_Tmp, 'DisplayName', 'Cooling Water Temp In [°C]')       % KVGT1 - Temp out of cooling HE
% plot(time, outCooler_Tmp, 'DisplayName', 'Cooling Water Temp Out [°C]')     % KVGT2 - Temp out of cooling HE
% plot(time, inHE_Tmp, 'DisplayName', 'Water Temp In HE [°C]')                % OTGT1 - Temp into main HE (water-glycol)
% plot(time, outHE_Tmp, 'DisplayName', 'Water Temp Out HE [°C]')              % OTGT2 - Temp out main HE (water-glycol) 
% plot(time, outHE_Tmp2, 'DisplayName', 'Water Temp Out HE 2 [°C]')           % OTGT3 - Temp out main HE - further down tube (water-glycol)  
% plot(time, mixValveCool_Tmp, 'DisplayName', 'Cool In Mixing Valve [°C]')    % OTGT4 - Cold Inlet Mixing Valve 
% plot(time, mixValveHot_Tmp, 'DisplayName', 'Hot In Mixing Valve [°C]')      % OTGT6 - Hot Inlet Mixing Valve 
plot(time, mixValve_Per, '-', 'Color', [0 0.4470 0.7410], 'DisplayName', 'Mixing Valve [%]', 'LineWidth', 2)         % OTSV1 - Mixing Valve Percent (ref.)
plot(time, outGearbox_Tmp, '-', 'Color', [0.9290 0.6940 0.1250], 'DisplayName', 'Oil Temp Out [°C]', 'LineWidth', 2)  % TV11 - Outlet Gearbox Temp 
plot(time, inGearbox_Tmp, '-', 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'Oil Temp In [°C]', 'LineWidth', 2)                % TV12 - Inlet Gearbox Temp
plot(time, inGearbox_Flw, '-', 'Color', [0.4660 0.6740 0.1880], 'DisplayName', 'Oil Flow [L/min.]', 'LineWidth', 2)                   % QV11 - Inlet Gearbox Flow
% plot(time, ref_Tmp, 'DisplayName', 'Reference Temp [°C]')                   % TV_Cmd - Reference TV12 Temp
% plot(time, mixValve_REAL, '-', 'Color', [0.3010 0.7450 0.9330], 'DisplayName', 'Real Mixing Valve [%]', 'LineWidth', 2)                % OTSV1_FB_LREAL - REAL Mixing Valve Percent
% plot(time, coolValve_Per, 'DisplayName', 'Oil Temp In [°C]')                % KVSV1 - Feedback for coolant valve  
yline(equil, '--k', 'DisplayName', 'Reference - 48 °C', 'LineWidth', 2)

ylabel('Temperature, Percentage, and Flow')
ylim([5 (round(max(outGearbox_Tmp), -1)+10)]);  
yticks(5:5:(round(max(outGearbox_Tmp), -1))+10);
ax.YColor = [0 0.4470 0.7410];

yyaxis right
% plot(time, torque, ':', 'Color',[0.4940 0.1840 0.5560], 'DisplayName', 'Torque [Nm]', 'LineWidth', 2)
plot(time, flow, ':', 'Color',[0.4940 0.1840 0.5560], 'DisplayName', 'Oil Pump RPM', 'LineWidth', 2)
ylabel('Oil Pump Speed [RPM]')
ylim([0 3500]);
yticks(0:350:3500);
ax.YColor = [0.4940 0.1840 0.5560];

title(['Flow Disturbance - ', num2str(flowRef), ' RPM and 500 Nm ', ctrl])
xlabel('Time [min.]')
xlim([0 max(time)])

legend()
hold off

%% FLOW - SCANIA

%% --- 350 RPM ---

T_data_all = readtable("H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\Scania\Scania_350.csv");
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(Manually Tuned)';
flowRef = 350;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(.43*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% --- 700 RPM ---

T_data_all = readtable("H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\Scania\Scania_700.csv");
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(Manually Tuned)';
flowRef = 700;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(.35*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% --- 1050 RPM ---

T_data_all = readtable("H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\Scania\Scania_1050.csv");
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(Manually Tuned)';
flowRef = 1050;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(.45*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% --- 1400 RPM ---

T_data_all = readtable("H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\Scania\Scania_1400.csv");
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(Manually Tuned)';
flowRef = 1400;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(.45*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% --- 2100 RPM ---

T_data_all = readtable("H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\Scania\Scania_2100.csv");
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(Manually Tuned)';
flowRef = 2100;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(.45*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% --- 2450 RPM ---

T_data_all = readtable("H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\Scania\Scania_2450.csv");
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(Manually Tuned)';
flowRef = 2450;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(.42*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% --- 2800 RPM ---

T_data_all = readtable("H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\Scania\Scania_2800.csv");
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(Manually Tuned)';
flowRef = 2800;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(.45*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% --- 3150 RPM ---

T_data_all = readtable("H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\Scania\Scania_3150.csv");
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(Manually Tuned)';
flowRef = 3150;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(.45*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% --- 3500 RPM ---

T_data_all = readtable("H:\Shared drives\Scania Thesis\01. Controller Design\Data\01 Test Data (T7 Apr.25)\Clean Data\Flow\Scania\Scania_3500.csv");
T_data_all.Properties.VariableNames = {'X1', 'KVTG1', 'X2', 'KVTG2', 'X3', 'OTGT1', 'X4', 'OTGT2', 'X5', 'OTGT3', 'X6', 'OTGT4', 'X7', 'OTGT6', 'X8', 'OTSV1', 'X9', 'QV11', 'X10', 'TV11', 'X11', 'TV12', 'X12', 'TV_Cmd', 'X13', 'OTSV1_FB_LREAL', 'X14', 'OTGP1', 'X15', 'KVSV1'};

Ts = .002;
ctrl = '(Manually Tuned)';
flowRef = 3500;
torqueRef = 500;
simStart = double(0*(60*(1/Ts)))+1;
torqueStart = double(.43*(60*(1/Ts))); % .13 min start
torqueEnd = double(torqueStart + (5*60*(1/Ts)));

%% ////////////// PLOT FOR FLOW: ARMAX ///////////////

% Raw data plotting variables (For all except 2100)
inCooler_Tmp = table2array(T_data_all(simStart:end, 2));           % KVGT1 - Temp out of cooling HE  
outCooler_Tmp = table2array(T_data_all(simStart:end, 4));          % KVGT2 - Temp out of cooling HE
inHE_Tmp = table2array(T_data_all(simStart:end, 6));               % OTGT1 - Temp into main HE (water-glycol)
outHE_Tmp = table2array(T_data_all(simStart:end, 8));              % OTGT2 - Temp out main HE (water-glycol) 
outHE_Tmp2 = table2array(T_data_all(simStart:end, 10));             % OTGT3 - Temp out main HE - further down tube (water-glycol)      
mixValveCool_Tmp = table2array(T_data_all(simStart:end, 12));       % OTGT4 - Cold Inlet Mixing Valve 
mixValveHot_Tmp = table2array(T_data_all(simStart:end, 14));        % OTGT6 - Hot Inlet Mixing Valve 
mixValve_Per = table2array(T_data_all(simStart:end, 16));           % OTSV1 - Mixing Valve Percent (ref.)
inGearbox_Flw = table2array(T_data_all(simStart:end, 18));          % QV11 - Inlet Gearbox Flow
outGearbox_Tmp = table2array(T_data_all(simStart:end, 20));        % TV11 - Outlet Gearbox Temp  
inGearbox_Tmp = table2array(T_data_all(simStart:end, 22));         % TV12 - Inlet Gearbox Temp
ref_Tmp = table2array(T_data_all(simStart:end, 24));               % TV_Cmd - Reference TV12 Temp
mixValve_REAL = table2array(T_data_all(simStart:end, 26));         % OTSV1_FB_LREAL - REAL Mixing Valve Percent
inHE_Pre = table2array(T_data_all(simStart:end, 28));              % OTGP1 - Pressure out of mixing valve
coolValve_Per = table2array(T_data_all(simStart:end, 30));         % KVSV1 - Coolant Valve Percent    

total_time = Ts*(length(mixValve_Per)-1);
time = (0:Ts:total_time)./60;
equil = 48;

% torque = zeros(1,size(time,2));
% for i=torqueStart:torqueEnd
%     torque(i) = torqueRef;
% end

flow = zeros(1,size(time,2))+1750;
for i=torqueStart:torqueEnd
    flow(i) = flowRef;
end

figure(1)
clf
hold on
ax = gca;
ax.FontSize = 18;
yyaxis left
% plot(time, inCooler_Tmp, 'DisplayName', 'Cooling Water Temp In [°C]')       % KVGT1 - Temp out of cooling HE
% plot(time, outCooler_Tmp, 'DisplayName', 'Cooling Water Temp Out [°C]')     % KVGT2 - Temp out of cooling HE
% plot(time, inHE_Tmp, 'DisplayName', 'Water Temp In HE [°C]')                % OTGT1 - Temp into main HE (water-glycol)
% plot(time, outHE_Tmp, 'DisplayName', 'Water Temp Out HE [°C]')              % OTGT2 - Temp out main HE (water-glycol) 
% plot(time, outHE_Tmp2, 'DisplayName', 'Water Temp Out HE 2 [°C]')           % OTGT3 - Temp out main HE - further down tube (water-glycol)  
% plot(time, mixValveCool_Tmp, 'DisplayName', 'Cool In Mixing Valve [°C]')    % OTGT4 - Cold Inlet Mixing Valve 
% plot(time, mixValveHot_Tmp, 'DisplayName', 'Hot In Mixing Valve [°C]')      % OTGT6 - Hot Inlet Mixing Valve 
plot(time, mixValve_Per, '-', 'Color', [0 0.4470 0.7410], 'DisplayName', 'Mixing Valve [%]', 'LineWidth', 2)         % OTSV1 - Mixing Valve Percent (ref.)
plot(time, outGearbox_Tmp, '-', 'Color', [0.9290 0.6940 0.1250], 'DisplayName', 'Oil Temp Out [°C]', 'LineWidth', 2)  % TV11 - Outlet Gearbox Temp 
plot(time, inGearbox_Tmp, '-', 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'Oil Temp In [°C]', 'LineWidth', 2)                % TV12 - Inlet Gearbox Temp
plot(time, inGearbox_Flw, '-', 'Color', [0.4660 0.6740 0.1880], 'DisplayName', 'Oil Flow [L/min.]', 'LineWidth', 2)                   % QV11 - Inlet Gearbox Flow
% plot(time, ref_Tmp, 'DisplayName', 'Reference Temp [°C]')                   % TV_Cmd - Reference TV12 Temp
% plot(time, mixValve_REAL, '-', 'Color', [0.3010 0.7450 0.9330], 'DisplayName', 'Real Mixing Valve [%]', 'LineWidth', 2)                % OTSV1_FB_LREAL - REAL Mixing Valve Percent
% plot(time, coolValve_Per, 'DisplayName', 'Oil Temp In [°C]')                % KVSV1 - Feedback for coolant valve  
yline(equil, '--k', 'DisplayName', 'Reference - 48 °C', 'LineWidth', 2)

ylabel('Temperature, Percentage, and Flow')
ylim([5 (round(max(outGearbox_Tmp), -1)+10)]);  
yticks(5:5:(round(max(outGearbox_Tmp), -1))+10);
ax.YColor = [0 0.4470 0.7410];

yyaxis right
% plot(time, torque, ':', 'Color',[0.4940 0.1840 0.5560], 'DisplayName', 'Torque [Nm]', 'LineWidth', 2)
plot(time, flow, ':', 'Color',[0.4940 0.1840 0.5560], 'DisplayName', 'Oil Pump RPM', 'LineWidth', 2)
ylabel('Oil Pump Speed [RPM]')
ylim([0 3850]);
yticks(0:350:3850);
ax.YColor = [0.4940 0.1840 0.5560];

title(['Flow Disturbance - ', num2str(flowRef), ' RPM and 500 Nm ', ctrl])
xlabel('Time [min.]')
xlim([0 max(time)])

legend()
hold off
