clc
%clear all

%% Data extraction

T = readtable('C:\Users\Blane Brinkley\OneDrive\Documents\Scania_GearboxOilTestBed\thesis_Data\testData_clean_011922.csv');
T.Properties.VariableNames = {'KVSV1', 'OTGT1', 'OTGT6', 'OTSV1', 'QV11', 'TV11', 'TV12', 'TV_Cmd'};
T = flip(T); % Data recorded backward

delay = 30; % 15 sec delay
Ts = 0.002;   % 2ms

train_low = 330000;
train_high = 720000;
test_low = 1200000;
test_high = 1560000;

desired_rate = .10; % [sec]
sample_rate = desired_rate/Ts;

mini_batch_train = T(train_low:sample_rate:train_high, :);
mini_batch_test  = T(test_low:sample_rate:test_high, :);

u_train = table2array(mini_batch_train(:, 4));   %OTSV1 (4) - 3-Way valve percentage
y_train = table2array(mini_batch_train(:, 7));   %TV12 (7) - Temp into gearbox

u_test = table2array(mini_batch_test(:, 4));   %OTSV1 (4) - 3-Way valve percentage
y_test = table2array(mini_batch_test(:, 7));   %TV12 (7) - Temp into gearbox

%max temp allowed is 100 degrees
%not pressurized system

% Hybrid Gearbox being used:
coolValve_percent = table2array(T(:, 1)); %KVSV1 (1) - Valve percentage for coolant temp
intoHE_temp = table2array(T(:, 2)); %OTGT1 (2)- Heat of fluid leading into heat exchanger
heater_temp = table2array(T(:, 3)); %OTGT6 (3) - Heater temp leading into 3-way valve
mainValve_percent = table2array(T(:, 4));%OTSV1 (4) - 3-Way valve percentage
intoGearbox_flow = table2array(T(:, 5)); %QV11 (5) - main valve percentage opening
intoGearbox_temp = table2array(T(:, 6)); %TV11 (6) - Temp from outlet of gearbox
outletGearbox_temp = table2array(T(:, 7)); %TV12 (7) - Temp into gearbox
refTemp = table2array(T(:, 8)); %TV_cmd (8) - reference temperature for closed loop


total_time = Ts*(length(coolValve_percent)-1);
time = (0:Ts:total_time)./60;


%% Plotting Initial Data

figure(1)
clf
hold on
plot(time, refTemp, 'DisplayName', 'Temp Reference (C)')
plot(time, mainValve_percent, 'DisplayName', 'Valve Percentage (0-100%)')
plot(time, intoHE_temp, 'DisplayName', 'Into HE Temp (C)')
plot(time, outletGearbox_temp, 'DisplayName', 'Gearbox Outlet Temp (C)')

plot(time, intoGearbox_flow, 'DisplayName', 'Gearbox Flow (L/min)')
%plot(time, coolValve_percent, 'DisplayName', 'Coolant Valve Percentage')
%plot(time, heater_temp, 'DisplayName', 'Heater Temp (C)')

%rectangle('position', [(train_low*Ts/60), 0, (train_high*Ts/60)-(train_low*Ts/60), 100], 'EdgeColor', 'b', 'LineStyle', '--')
%rectangle('position', [(test_low*Ts/60), 0, (test_high*Ts/60)-(test_low*Ts/60), 100], 'EdgeColor', 'r', 'LineStyle', '--')
title('Closed Loop Controller Plotting')
xlabel('Time (minutes)')
legend()
hold off

%% Training Data
data_train = iddata(y_train, u_train, desired_rate);

data_train.Name = 'mini_batch_train';
data_train.TimeUnit = 'seconds';

data_train.InputName = 'Valve %';   data_train.InputUnit = 'Percentage';
data_train.OutputName = 'Gearbox Input Temp';   data_train.OutputUnit = 'Celsius';

% Filter Specs
data_train_filter = idfilt(data_train, [0 1]); %butter worth bandpass filter with 0-1 hz 
data_train_filter = data_train_filter(25/desired_rate:end, :, :) 

%% Testing Data
data_test = iddata(y_test, u_test, desired_rate);

data_test.Name = 'mini_batch_test';
data_test.TimeUnit = 'seconds';

data_test.InputName = 'Valve %';   data_test.InputUnit = 'Percentage';
data_test.OutputName = 'Gearbox Input Temp';   data_test.OutputUnit = 'Celsius';

data_test_filter = idfilt(data_test, [0 1]);
data_test_filter = data_test_filter(25/desired_rate:end, :, :) 

%% Print training data
get(data_train_filter)
advice(data_train_filter, 'nonlinearity')

figure(2)
clf
hold on
plot(data_train, 'r', data_train_filter, 'b')
%plot(data_train_filter, 'b')
hold off

figure(3)
clf
hold on
plot(data_test, 'r')
plot(data_test, 'r', data_test_filter, 'b')
hold off

%% Functions For Running Different models 

%sysTF = tfFunction(data_train_filter, data_test_filter, delay)
% procEst = processFunction(data_train_filter, data_test_filter, delay)
% [ssEst, ssTFEst] = ssFunction(data_train_filter, data_test_filter, delay)
% arxEst = arxFunction(data_train_filter, data_test_filter, delay)
% armaxEst = armaxFunction(data_train_filter, data_test_filter, delay)
% bjEst = bjFunction(data_train_filter, data_test_filter, delay)

%% Estimate TF Model

function sysTF = tfFunction(data_train_filter, data_test_filter, delay)
    opt = tfestOptions;
    opt.InitializeMethod = 'all';
    opt.Focus = 'prediction';
    opt.SearchOptions.MaxIterations = 100; 
    opt.Display = 'on';

    np = 4;          % Num of pole
    nz = 3;          % Num of zero
    iodelay = delay;   % In/Out delay

    sysTF = tfest(data_train_filter, np, nz, iodelay, opt)

    %advice(sysTF, data_test_filter)

    figure(4)
    compare(data_test_filter, sysTF)
    
    figure(5)
    resid(data_test_filter,sysTF)
end 

%% Process Model Estimation with Disturbance Model

function procEst = processFunction(data_train_filter, data_test_filter, delay)
    % sysInit = idproc('P3D','TimeUnit','seconds');
    % 
    % sysInit.Structure.Kp.Value    = 1;
    % sysInit.Structure.Kp.Minimum  = 0;
    % sysInit.Structure.Tp1.Value   = 1;
    % sysInit.Structure.Tp1.Maximum = 10;
    % sysInit.Structure.Td.Value    = 15;
    % sysInit.Structure.Td.Minimum  = 10;
    % sysInit.Structure.Td.Maximum  = 20;
    % 
    % opt = procestOptions('DisturbanceModel','ARMA1');
    %procEst = procest(data_train,sysInit,opt);
    procEst = procest(data_train_filter, 'P3DU');

    figure(6)
    compare(data_test_filter, procEst)
    
    figure(7)
    resid(data_test_filter, procEst)
end

%% State Space Model

function [ssEst, ssTFEst] = ssFunction(data_train_filter, data_test_filter, delay)
    opt = ssestOptions;
    opt.InitializeMethod = 'auto';
    opt.Focus = 'prediction';
    opt.SearchOptions.MaxIterations = 50; 
    opt.Display = 'on';

    nx = 4
    ssEst = ssest(data_train_filter, nx, opt)

    [b, a] = ss2tf(sysSS.A, sysSS.B, sysSS.C, sysSS.D);
    ssTFEst = tf(b, a)
    
    figure(8)
    compare(data_test_filter, ssESt, ssTFEst)
    
    figure(9)
    resid(data_test_filter,ssEst)
    
    figure(10)
    resid(data_test_filer, ssTFest)

end
%% ARX

function arxEst = arxFunction(data_train_filter, data_test_filter, delay)

    % na = order of A(q) polynomial
    % na = order of B(q) polynomial
    % nk = input-output delay
    arxEst = arx(data_train_filter, [1 2 10]) %[na nb nk] number of coefficients for the polynomial
   
    figure(11)
    compare(data_test_filter, ssESt, ssTFEst)
    
    figure(12)
    resid(data_test_filter,ssEst)
    
end

%% ARMAX (includes a C polynomial for noise)

function armaxEst = armaxFunction(data_train_filter, data_test_filter, delay)
    
    armaxEst = armax(data_train_filter, [2 2 1 1]) %[na nb nc nk]
    
    figure(13)
    compare(data_test_filter, ssESt, ssTFEst)
    
    figure(14)
    resid(data_test_filter,ssEst)
end
%% Box-jenkins (BJ)

function bjEst = bjFunction(data_train_filter, data_test_filter, delay)
    bjEst = bj(data_train_filter, [2 2 2 2 1]) %[nb nc nd nf nk]

    figure(15)
    compare(data_test_filter, ssESt, ssTFEst)
    
    figure(16)
    resid(data_test_filter,ssEst)
end







