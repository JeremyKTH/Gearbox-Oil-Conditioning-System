clc
%clear all

%% DATA EXTRACTION 

%======JEREMY DATA=======
% SYSTEM ID DATA: 40 - 70 - 40 - 10 - 40  ( 30%) - D_47414
%T_data = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\2022-Feb-24\Partition Data\D_47414.csv');
%T_data.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11'};

%======BLANE (ALL) DATA===========
% SYSTEM ID DATA: 40 - 70 - 40 - 10 - 40  ( 30%) - D_47414
T_data = readtable('E:\Users\Blane Brinkley\OneDrive\Documents\Scania_GearboxOilTestBed\thesis_Data\systemID_testing_2.24\segmented_Data\D_47414_all.csv');
T_data.Properties.VariableNames = {'OTSV1', 'QV11', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'KVSV1', 'KVSV1_FB', 'KVGT2', 'TV12'};

delay = 20; % 15 sec delay
Ts = 0.2;   % 200 ms --> base time from when Jeremy converted

% Raw data plotting variables
mainValve_Per = table2array(T_data(:, 1));          % OTSV1 (1) - Mixing valve percentage
intoGearbox_Flow = table2array(T_data(:, 2));       % QV11 (2) - Inlet Gearbox Flow                                           
outGearbox_Temp = table2array(T_data(:, 3));        % TV11 (3) - Outlet Gearbox Temp                                           
intoHE_Temp = table2array(T_data(:, 4));            % OTGT1 (4) - Inlet HE Temp                                          
outHE_Temp = table2array(T_data(:, 5));             % OTGT2 (5) - Outlet HE Temp                                          
intoCooler_Temp = table2array(T_data(:, 6));        % OTGT3 (6) - Inlet Cooler Temp                                           
intoMix_cold_Temp = table2array(T_data(:, 7));      % OTGT4 (7) - Cold Inlet Mixing Valve                                           
intoMix_hot_Temp = table2array(T_data(:, 8));       % OTGT6 (8) - Heat Inlet Mixing Valve                                          
coolantValve_Per = table2array(T_data(:, 9));       % KVSV1 (9) - Coolant Valve Percentage                                          
coolantValve_FB = table2array(T_data(:, 10));       % KVSV1_FB (10) - Feedback for coolant valve                                       
coolantSource_Temp = table2array(T_data(:, 11));    % KVGT2 (11) - Source temp for cooler                                          
intoGearbox_temp = table2array(T_data(:, 12));      % TV12 (12) - Inlet Gearbox Temp

%=====Training data sets=====
T_test = flip(T_data); % Flipping the data for 

u_train = table2array(T_data(:, 1));   %OTSV1 (1) - 3-Way valve percentage
y_train = table2array(T_data(:, 12));   %TV12 (12) - Temp into gearbox

u_test = table2array(T_test(:, 1));   %OTSV1 (1) - 3-Way valve percentage
y_test = table2array(T_test(:, 12));   %TV12 (12) - Temp into gearbox

total_time = Ts*(length(u_train)-1); % (sec)
time = (0:Ts:total_time)./60;
%% Plotting Initial Data

figure(1)
clf
hold on
plot(time, mainValve_Per, 'DisplayName', 'Mixing valve percentage')
%plot(time, intoGearbox_Flow, 'DisplayName', 'Inlet Gearbox Flow)')
plot(time, outGearbox_Temp, 'DisplayName', 'Outlet Gearbox Temp')
plot(time, intoHE_Temp, 'DisplayName', 'Inlet HE Temp')
%plot(time, outHE_Temp, 'DisplayName', 'Outlet HE Temp')
%plot(time, intoCooler_Temp, 'DisplayName', 'Inlet Cooler Temp')
%plot(time, intoMix_cold_Temp, 'DisplayName', 'Cold Inlet Mixing Valve')
plot(time, intoMix_hot_Temp, 'DisplayName', 'Heat Inlet Mixing Valve ')
plot(time, coolantValve_Per, 'DisplayName', 'Coolant Valve Percentage')
%plot(time, coolantValve_FB, 'DisplayName', 'Feedback for coolant valve')
%plot(time, coolantSource_Temp, 'DisplayName', 'Source temp for cooler')
plot(time, intoGearbox_temp, 'DisplayName', 'Inlet Gearbox Temp')
title('Raw Data Plot')
xlabel('Time (minutes)')
legend()
hold off

%% TRAINING DATA
data_train = iddata(y_train, u_train, Ts);

data_train.Name = 'Data_Train';
data_train.TimeUnit = 'seconds';

data_train.InputName = 'Valve% (OTSV1)';   data_train.InputUnit = 'Percentage';
data_train.OutputName = 'Gearbox Input Temp (TV12)';   data_train.OutputUnit = 'Celsius';

data_train = detrend(data_train);

%=== FILTERING ===
fs = 1/Ts; % sample rate (Hz) 5 Hz = .2 sec
fn = fs/2;
fhc = .2/fn; % high cut-off frequency (Hz)
flc = 0/fn; % low cut-off frequency (Hz)
wn = [flc fhc];
N = 5; % filer order

%--- Original (idfilt) butterworth --
%data_train_filt_old = idfilt(data_train, [0 1]); %butter worth bandpass filter with 0-1 hz 
%data_train_filt.InputData = data_train.InputData;
%data_train_filter = data_train_filter(25/desired_rate:end, :, :) %Trimming

%--- Manual butterwroth ---

[z,p,k] = butter(N, fhc,'low'); 
%returns: zeros, poles, gain
%butter(order, bandpass, type)

sos = zp2sos(z,p,k);
figure(2)
phasedelay(sos,1024);
figure(3)
grpdelay(sos,1024)
[b,a] = zp2tf(z,p,k);

data_train_filt.OutputData = filtfilt(b, a, data_train.OutputData);

%--- Manual IIR Low Pass ---
%Niir = 5; %orderd of the filter (butter nom. = 5)
%iir = designfilt('lowpassiir','FilterOrder',Niir,'HalfPowerFrequency',fc,'SampleRate',fs) %Low-pass IIR filter
%data_train_filt.OutputData = filtfilt(iir, data_train.OutputData);

figure(4)
clf
hold on
plot(data_train, 'b', data_train_filt, 'r')
%plot(data_train_filter, 'b')
hold off

%% TESTING DATA
data_test = iddata(y_test, u_test, Ts);

data_test.Name = 'Data_Test';
data_test.TimeUnit = 'seconds';

data_test.InputName = 'Valve% (OTSV1)';   data_test.InputUnit = 'Percentage';
data_test.OutputName = 'Gearbox Input Temp (TV12)';   data_test.OutputUnit = 'Celsius';

data_test = detrend(data_test);

% Filter Specs
data_test_filt = idfilt(data_test, [0 1]);
data_test_filt.InputData = data_test.InputData;
%data_test_filter = data_test_filter(25/desired_rate:end, :, :) %Removing


figure(3)
clf
hold on
%plot(data_test, 'b')
plot(data_test, 'b', data_test_filt, 'r')
hold off

%% Functions For Running Different models 

sysTF = tfFunction(data_train_filt, data_test_filt, delay)
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

    np = 2;          % Num of pole
    nz = 1;          % Num of zero
    iodelay = delay;   % In/Out delay

    sysTF = tfest(data_train_filter, np, nz, iodelay, opt)
    
    advice(sysTF, data_test_filter)
    sysTF.Report.Fit
    
    % generate reference output
        y_tf_ref = data_train.y;
    % generate estimated output
        y_tf_est = sim(sysTF, data_train(:, [], :));
    % select cost function
        cost_func = 'NRMSE';

    % Goodness of Fit
    fit = goodnessOfFit(y_tf_est.y, y_tf_ref, cost_func)
    
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







