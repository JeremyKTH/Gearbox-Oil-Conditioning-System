clc
clear all

%% DATA EXTRACTION 

%======JEREMY DATA=======
% SYSTEM ID DATA: 40 - 70 - 40 - 10 - 40  ( 30%) - D_47414
%T_data = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\2022-Feb-24\Partition Data\D_47414.csv');
%T_data.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11'};

%======BLANE (ALL) DATA===========
% SYSTEM ID DATA: 40 - 70 - 40 - 10 - 40  ( 30%) - D_47414
T_data_all = readtable('E:\Users\Blane Brinkley\OneDrive\Documents\Scania_GearboxOilTestBed\thesis_Data\systemID_testing_2.24\segmented_Data\D_49404_all.csv');
T_data_all.Properties.VariableNames = {'OTSV1', 'QV11', 'TV11', 'OTGT1', 'OTGT2', 'OTGT3', 'OTGT4', 'OTGT6', 'KVSV1', 'KVSV1_FB', 'KVGT2', 'TV12'};

delay = 20; % 15 sec delay
Ts_all = 0.4;   % 200 ms --> base time from when Jeremy converted

% Raw data plotting variables
mainValve_Per = table2array(T_data_all(:, 1));          % OTSV1 (1) - Mixing valve percentage
intoGearbox_Flow = table2array(T_data_all(:, 2));       % QV11 (2) - Inlet Gearbox Flow                                           
outGearbox_Temp = table2array(T_data_all(:, 3));        % TV11 (3) - Outlet Gearbox Temp                                           
intoHE_Temp = table2array(T_data_all(:, 4));            % OTGT1 (4) - Inlet HE Temp                                          
outHE_Temp = table2array(T_data_all(:, 5));             % OTGT2 (5) - Outlet HE Temp                                          
intoCooler_Temp = table2array(T_data_all(:, 6));        % OTGT3 (6) - Inlet Cooler Temp                                           
intoMix_cold_Temp = table2array(T_data_all(:, 7));      % OTGT4 (7) - Cold Inlet Mixing Valve                                           
intoMix_hot_Temp = table2array(T_data_all(:, 8));       % OTGT6 (8) - Heat Inlet Mixing Valve                                          
coolantValve_Per = table2array(T_data_all(:, 9));       % KVSV1 (9) - Coolant Valve Percentage                                          
coolantValve_FB = table2array(T_data_all(:, 10));       % KVSV1_FB (10) - Feedback for coolant valve                                       
coolantSource_Temp = table2array(T_data_all(:, 11));    % KVGT2 (11) - Source temp for cooler                                          
intoGearbox_temp = table2array(T_data_all(:, 12));      % TV12 (12) - Inlet Gearbox Temp

%% Training/Testing data sets
%T_data = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\Ts = 0.4\Training_Data\D_48404_train.csv');
%T_train.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1'};
T_train = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\Ts = 1\Training_Data\D_49404_train.csv');
T_train.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11'};

%T_data = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\Ts = 0.4\Training_Data\D_48404_test.csv');
%T_test.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11', 'OTGT1'};
T_test = readtable('H:\Shared drives\Scania Thesis\Code\Test Data\Ts = 1\Testing_Data\D_49404_test.csv');
T_test.Properties.VariableNames = {'OTSV1', 'TV12', 'TV11'};

%Ts = 0.4;   % 200 ms --> base time from when Jeremy converted
Ts = 1;

u_train = table2array(T_train(:, 1));   %OTSV1 (1) - 3-Way valve percentage
y_train = table2array(T_train(:, 2));   %TV12 (12) - Temp into gearbox

u_test = table2array(T_test(:, 1));   %OTSV1 (1) - 3-Way valve percentage
y_test = table2array(T_test(:, 2));   %TV12 (12) - Temp into gearbox

total_time = Ts_all*(length(mainValve_Per)-1); % (sec)
time = (0:Ts_all:total_time)./60;
%% Plotting Initial Data

figure(1)
clf
hold on
plot(time, mainValve_Per, 'DisplayName', 'Mixing valve percentage')
plot(time, intoGearbox_Flow, 'DisplayName', 'Inlet Gearbox Flow)')
plot(time, outGearbox_Temp, 'DisplayName', 'Outlet Gearbox Temp')
plot(time, intoHE_Temp, 'DisplayName', 'Inlet HE Temp')
%plot(time, outHE_Temp, 'DisplayName', 'Outlet HE Temp')
%plot(time, intoCooler_Temp, 'DisplayName', 'Inlet Cooler Temp')
plot(time, intoMix_cold_Temp, 'DisplayName', 'Cold Inlet Mixing Valve')
plot(time, intoMix_hot_Temp, 'DisplayName', 'Heat Inlet Mixing Valve ')
plot(time, coolantValve_Per, 'DisplayName', 'Coolant Valve Percentage')
%plot(time, coolantValve_FB, 'DisplayName', 'Feedback for coolant valve')
%plot(time, coolantSource_Temp, 'DisplayName', 'Source temp for cooler')
plot(time, intoGearbox_temp, 'DisplayName', 'Inlet Gearbox Temp')
title('Raw Data Plot (+/- 50%)')
xlabel('Time (minutes)')
legend()
hold off

%% TRAINING DATA
data_train = iddata(y_train, u_train, Ts);

data_train.Name = 'Data_Train';
data_train.TimeUnit = 'seconds';
data_train.InputName = 'Valve% (OTSV1)';   data_train.InputUnit = 'Percentage';
data_train.OutputName = 'TV12';   data_train.OutputUnit = 'Celsius';

%data_train = detrend(data_train);

%=== FILTERING ===
fs = 1/Ts; % sample rate (Hz) 5 Hz = .2 sec
fn = fs/2;
fhc = .2/fn; % high cut-off frequency (Hz)
flc = 0/fn; % low cut-off frequency (Hz)
wn = [flc fhc];
N = 5; % filer order

%--- Original (idfilt) butterworth --
data_train_filt = idfilt(data_train, [0 1]); %butter worth bandpass filter with 0-1 hz 
data_train_filt.InputData = data_train.InputData;
%data_train_filter = data_train_filter(25/desired_rate:end, :, :) %Trimming

%--- Manual butterwroth ---
[z,p,k] = butter(N, fhc,'low'); 
%returns: zeros, poles, gain
%butter(order, bandpass, type)

sos = zp2sos(z,p,k);
%figure(2)
%phasedelay(sos,1024);
%figure(3)
%grpdelay(sos,1024)
[b,a] = zp2tf(z,p,k);

data_train_filt.OutputData = filtfilt(b, a, data_train.OutputData);

%--- Manual IIR Low Pass ---
%Niir = 5; %orderd of the filter (butter nom. = 5)
%iir = designfilt('lowpassiir','FilterOrder',Niir,'HalfPowerFrequency',fc,'SampleRate',fs) %Low-pass IIR filter
%data_train_filt.OutputData = filtfilt(iir, data_train.OutputData);

% figure(4)
% clf
% hold on
% plot(data_train, 'b', data_train_filt, 'r')
% hold off

%% TESTING DATA
data_test = iddata(y_test, u_test, Ts);

data_test.Name = 'Data_Test';
data_test.TimeUnit = 'seconds';
data_test.InputName = 'Valve% (OTSV1)';   data_test.InputUnit = 'Percentage';
data_test.OutputName = 'TV12';   data_test.OutputUnit = 'Celsius';

%data_test = detrend(data_test);

% Filter Specs
data_test_filt = idfilt(data_test, [0 1]);
data_test_filt.InputData = data_test.InputData;
%data_test_filter = data_test_filter(25/desired_rate:end, :, :) %Removing

data_test_filt.OutputData = filtfilt(b, a, data_test.OutputData);

% figure(5)
% clf
% hold on
% plot(data_test, 'b', data_test_filt, 'r')
% hold off

%% Detrending Data
T_train = getTrend(data_train);
T_test = getTrend(data_test);
T_train.InputOffset = 40;
T_train.OutputOffset = 48;
T_test.InputOffset = 40;
T_test.OutputOffset = 48;

data_train = detrend(data_train);
data_test = detrend(data_test);

cost_func = 'NRMSE';

%% Functions For Running Different models 

sysTF = tfFunction(data_train_filt, data_test_filt, data_test, data_train, delay)
% procEst = processFunction(data_train_filt, data_test_filt,data_test, data_train, delay)
% [sysSS, ssTFEst] = ssFunction(data_train_filt, data_test_filt,data_test, data_train, delay)
%arxEst = arxFunction(data_train_filt, data_test_filt, data_test, data_train, delay)
%sysARMAX = armaxFunction(data_train_filt, data_test_filt,data_test, data_train, delay);
%bjEst = bjFunction(data_train_filt, data_test_filt,data_test, data_train, delay)
    
%% Random NRMSE Generation
%     nrmsePlot =    [1, 0.2144;
%                     2, 0.2267;  
%                     3, 0.2135; 
%                     4, 0.2132];
%     figure(6)
%     plot(nrmsePlot(:,1), nrmsePlot(:,2))
%     title('SS Model Selection')
%     xlabel('Fitting Order')
%     ylabel('NRMSE')
%     xticks(0:1:4)    
        
%% Estimate Transfer Function
function sysTF = tfFunction(data_train_filter, data_test_filter, data_test, data_train, delay)
    opt = tfestOptions;
    opt.InitializeMethod = 'all';
    opt.Focus = 'prediction';
    opt.SearchOptions.MaxIterations = 1000; 
    opt.Display = 'on';
    %opt.DistrubanceModel = 'estimate'

    np = 1;          % Num of pole
    nz = 1;          % Num of zero
    iodelay = delay;   % In/Out delay

    sysTF = tfest(data_train, np, nz, iodelay, opt)
    
    sysTF.Report.Fit
    advice(sysTF, data_test)
    
    % generate reference output
    y_tf_ref = data_test.y;
    % generate estimated output
    y_tf_est = sim(sysTF, data_test(:, [], :));
    % select cost function
    cost_func = 'NRMSE';

    % Goodness of Fit
    fit = goodnessOfFit(y_tf_est.y, y_tf_ref, cost_func)
    
    opt = compareOptions('InitialCondition','e');   
%     figure (6)
%     compare(data_train_filter, sysTF, opt);
    
    figure(7)
    compare(data_test, sysTF, opt);
    
    figure(8)
    resid(data_test,sysTF);
    
%     nrmsePlot =    [1, .3758;
%                     2, .379;  
%                     3, .2583; 
%                     4, .3868];
%     figure(9)
%     plot(nrmsePlot(:,1), nrmsePlot(:,2))
%     title('TF Model Selection')
%     xlabel('Fitting Order')
%     ylabel('NRMSE')
%     xticks(0:1:4)
end 


%% State Space Model

function [sysSS, ssTFEst] = ssFunction(data_train_filter, data_test_filter, data_test, data_train, delay)
%     opt = ssestOptions;
%     opt.InitializeMethod = 'auto';
%     opt.Focus = 'prediction';
%     opt.SearchOptions.MaxIterations = 50; 
%     opt.Display = 'on';

    CoarseOrderSearch = [1, 2, 3, 4, 5, 6, 7, 8];
    %sysSS = ssest(data_train, CoarseOrderSearch, 'Ts', data_train.Ts)
    %sysSS = ssregest(data_train, CoarseOrderSearch, 'Ts', data_train.Ts)
    sysSS = n4sid(data_train, CoarseOrderSearch)
    
    %sysSS = ssest(data_train_filter, 1, 'Ts' , Ts)
    
    advice(sysSS, data_test)
    
    % generate reference output
    y_ss_ref = data_test.y;
    % generate estimated output
    y_ss_est = sim(sysSS, data_test(:, [], :));
    % select cost function
    cost_func = 'NRMSE';
    % Goodness of Fit
    fit = goodnessOfFit(y_ss_est.y, y_ss_ref, cost_func)
    
     [b, a] = ss2tf(sysSS.A, sysSS.B, sysSS.C, sysSS.D);
     ssTFEst = tf(b, a)
    
    % Specify an initial condition of zero to match the initial condition that goodnessOfFit assumes
    opt = compareOptions('InitialCondition','e');   % 'z' : zero initial conditions
    
%     figure(6)     
%     compare(data_test_filter, ssTFEst, opt);

    
    figure(7)
    compare(data_test, sysSS, opt);

    figure(8)
    resid(data_test, sysSS);
    
%     nrmsePlot =    [1, 1;
%                     2, 1;  
%                     3, .424; 
%                     4, .4433];
%     figure(9)
%     plot(nrmsePlot(:,1), nrmsePlot(:,2))
%     title('SS Model Selection')
%     xlabel('Fitting Order')
%     ylabel('NRMSE')
%     xticks(0:1:4)
end
%% ARX

function arxEst = arxFunction(data_train_filter, data_test_filter, data_test, data_train, delay)

    % na = order of A(q) polynomial
    % na = order of B(q) polynomial
    % nk = input-output delay
    arxEst = arx(data_train, [3 1 1]) %[na nb nk] number of coefficients for the polynomial
   
    arxEst.Report.Fit
    advice(arxEst, data_test)
    
    % generate reference output
    y_tf_ref = data_test.y;
    % generate estimated output
    y_tf_est = sim(arxEst, data_test(:, [], :));
    % select cost function
    cost_func = 'NRMSE';

    % Goodness of Fit
    fit = goodnessOfFit(y_tf_est.y, y_tf_ref, cost_func)
    
    opt = compareOptions('InitialCondition','e');   
%     figure (6)
%     compare(data_train_filter, arxEst, opt);
    
    figure(7)
    compare(data_test, arxEst, opt);
    
    figure(8)
    resid(data_test, arxEst);
    
%     nrmsePlot =    [1, .2153;
%                     2, .2148;  
%                     3, .2136; 
%                     4, .2124];
%     figure(9)
%     plot(nrmsePlot(:,1), nrmsePlot(:,2))
%     title('ARX Model Selection')
%     xlabel('Fitting Order')
%     ylabel('NRMSE')
%     xticks(0:1:4)
end

%% ARMAX (includes a C polynomial for noise)

function sysARMAX = armaxFunction(data_train_filter, data_test_filter, data_test, data_train, delay)
    
    opt = armaxOptions;
    opt.Focus = 'prediction';
    opt.SearchOptions.MaxIterations = 1000;
    opt.SearchOptions.Tolerance = 1e-5;

    na = 2;
    nb = 1;
    nc = 1;
    nk = 1;
    
    sysARMAX = armax(data_train, [na nb nc nk], opt)

    advice(sysARMAX, data_test)
    
    % generate reference output
    y_armax_ref = data_test.y;
    % generate estimated output
    y_armax_est = sim(sysARMAX, data_test(:, [], :));
    % select cost function
    cost_func = 'NRMSE';
    % Goodness of Fit
    fit = goodnessOfFit(y_armax_est.y, y_armax_ref, cost_func)
    
    opt = compareOptions('InitialCondition','e');   
    figure(7)
    compare(data_test, sysARMAX, opt)
    
    figure(8)
    resid(data_test,sysARMAX)
    
    figure(8)
    resid(data_test_filter,arxEst);
    
%     nrmsePlot =    [1, .2117;
%                     2, .305;  
%                     3, .2145; 
%                     4, .2432];
%     figure(9)
%     plot(nrmsePlot(:,1), nrmsePlot(:,2))
%     title('ARMAX Model Selection')
%     xlabel('Fitting Order')
%     ylabel('NRMSE')
%     xticks(0:1:4)
end
%% Box-jenkins (BJ)

function bjEst = bjFunction(data_train_filter, data_test_filter, data_test, data_train, delay)
    opt = bjOptions;
    opt.Focus = 'prediction';
    opt.SearchOptions.MaxIterations = 1000;
    opt.SearchOptions.Tolerance = 1e-5;

    nb = 2; %num
    nc = 1; %num (noise)
    nd = 1; %den (noise)
    nf = 2; %den
    nk = 1;

    bjEst = bj(data_train, [nb nc nd nf nk], opt, 'IODelay', 20) %[nb nc nd nf nk]
    
    advice(bjEst, data_test)
    
    % generate reference output
    y_bj_ref = data_test.y;
    % generate estimated output
    y_bj_est = sim(bjEst, data_test(:, [], :));
    % select cost function
    cost_func = 'NRMSE';

    % Goodness of Fit
    fit = goodnessOfFit(y_bj_est.y, y_bj_ref, cost_func)
    
    opt = compareOptions('InitialCondition','e'); 
    figure(7)
    compare(data_test, bjEst, opt);
    
    figure(8)
    resid(data_test,bjEst)
    
%     nrmsePlot =    [1, .3868;
%                     2, .3892;  
%                     3, .3777; 
%                     4, .3803];
%     figure(9)
%     plot(nrmsePlot(:,1), nrmsePlot(:,2))
%     title('BJ Model Selection')
%     xlabel('Fitting Order')
%     ylabel('NRMSE')
%     xticks(0:1:4)
end

%% Process Model Estimation with Disturbance Model

function procEst = processFunction(data_train_filter, data_test_filter, data_test, data_train, delay)
    sysInit = idproc( 'P2Z','TimeUnit','seconds');
    %sysInit = idproc( 'P1D','TimeUnit','seconds', 'inputDelay', 20)
    % 
    % sysInit.Structure.Kp.Value    = 1;
    % sysInit.Structure.Kp.Minimum  = 0;
    % sysInit.Structure.Tp1.Value   = 1;
    % sysInit.Structure.Tp1.Maximum = 10;
    % sysInit.Structure.Td.Value    = 15;
    % sysInit.Structure.Td.Minimum  = 10;
    % sysInit.Structure.Td.Maximum  = 20;
    % 
    opt = procestOptions;
    opt.DisturbanceModel = 'Estimate'; %recommends 'Estimate' (old: 'ARMA1')
    %opt.Focus = 'simulation'; %Alternative = 'simulation', 'prediction' = one-step ahead
    %opt.InitialCondition = 'estimate'; %'zero' is initial condition
    
    %procEst = procest(data_train,sysInit,opt);
    procEst = procest(data_train, sysInit, opt);

    advice(procEst, data_test)
    
    % generate reference output
    y_bj_ref = data_test.y;
    % generate estimated output
    y_bj_est = sim(procEst, data_test(:, [], :));
    % select cost function
    cost_func = 'NRMSE';
    % Goodness of Fit
    fit = goodnessOfFit(y_bj_est.y, y_bj_ref, cost_func)
    
    opt = compareOptions('InitialCondition','e'); 
    figure(7)
    compare(data_test, procEst, opt);
    
    figure(8)
    resid(data_test,procEst)
    
%     nrmsePlot =    [1, .2105;
%                     2, .2077;  
%                     3, .2129]; 
%                     
%     figure(9)
%     plot(nrmsePlot(:,1), nrmsePlot(:,2))
%     title('Process Model Selection')
%     xlabel('Fitting Order')
%     ylabel('NRMSE')
%     xticks(0:1:4)
end






