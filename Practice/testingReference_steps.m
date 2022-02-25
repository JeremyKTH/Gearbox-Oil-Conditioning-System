clc
clear all

%Initial Notes...
%   ~12 min to heat up 

%% Desired Values (INPUT HERE)
equil = 70; % equilibrium point for valve openeing [%]
heatTime = 12; %time to heat-up [min]
coolDownTime = 20; % time to cooldown [min]
sampleRate = .1; % desired rate for sending reference [sec]
startTime = 20; % when to start the square wave[min]
period = 10; % period of square wave [min]

%% Conversions

heatTime = heatTime*60/sampleRate; %sec
startTime = startTime*60/sampleRate;
period = period*60/sampleRate;
coolDownTime = coolDownTime*60/sampleRate;

if equil-0 > 50
    counterMax = (100-equil)/10;
else
    counterMax = (equil-0)/10;
end

totalTime = (startTime+coolDownTime+(period*counterMax)+period)*sampleRate; %[sec]

time = 0:sampleRate:totalTime;

%% Step-up Reference Generation

% Heat up/cooldown
for i=1:size(time,2)
    if i <= heatTime
        valve(i) = equil;
    elseif i > (startTime + counterMax*period + period)
        valve(i) = 30;
    else
        valve(i) = equil;
    end
end

for n = 0:counterMax-1
    for i = (startTime + period*n)+1:(startTime + period*(n+1)+1)
        valve(i) = equil + (10*(n+1));
    end
end

figure(1)
clf 
plot(time./60, valve)
%plot(time, valve)
ylim([0 110])
xlim([0 totalTime/60])
%xlim([0 totalTime])
yline(equil, '--r')
title('Reference Input for Linearziation Validity (Step-Up)')
xlabel('Time [min]')
ylabel('Reference Temperature [°C]')

% m = [time; valve]';
% writematrix(m, 'systemID_ReferenceInput.csv');

%% Step-down Reference Generation

counterMax = (equil-30)/10
% Heat up/cooldown
for i=1:size(time,2)
    if i <= heatTime
        valve(i) = equil;
    elseif i > (startTime + counterMax*period + period)
        valve(i) = 30;
    else
        valve(i) = equil;
    end
end

for n = 0:counterMax-1
    for i = (startTime + period*n)+1:(startTime + period*(n+1)+1)
        valve(i) = equil - (10*(n+1));
    end
end

figure(2)
clf 
plot(time./60, valve)
%plot(time, valve)
ylim([0 110])
xlim([0 totalTime/60])
%xlim([0 totalTime])
yline(equil, '--r')
title('Reference Input for Linearziation Validity (Step-Down)')
xlabel('Time [min]')
ylabel('Reference Temperature [°C]')
% m = [time; valve]';
% writematrix(m, 'systemID_ReferenceInput.csv');

