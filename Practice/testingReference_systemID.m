clc
clear all

%Initial Notes...
%   ~12 min to heat up 

%% Desired Values (INPUT HERE)
equil = 55; % equilibrium point for valve openeing [%]
heatTime = 20; %time to heat-up [min]
coolDownTime = 20; % time to cooldown [min]
sampleRate = .1; % desired rate for sending reference [sec]
startTime = 40; % when to start the square wave[min]
period = 80; % period of square wave [min]

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

totalTime = (startTime+coolDownTime+(period*counterMax))*sampleRate; %[sec]

time = 0:sampleRate:totalTime;

%% Generating reference 

% Heat up/cooldown
for i=1:size(time,2)
    if i <= heatTime
        valve(i) = 100;
    elseif i >= (startTime + counterMax*period + 1)
        valve(i) = 30;
    else
        valve(i) = equil;
    end
end

for n = 0:counterMax-1
    for i = (startTime + period*n)+1:(startTime + period*(n+1)+1)
        if i <= startTime + period*(n + 1/4)
            valve(i) = equil + (10*(n+1));
        elseif i > (startTime + period*(n+1/2)) && i<= (startTime + period*(n+3/4))
            valve(i) = equil-(10*(n+1));
        end
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
title('Reference Input for System Identification')
xlabel('Time [min]')
ylabel('Input - Valve Percentage [%]')

m = [time; valve]';
writematrix(m, 'systemID_ReferenceInput.csv');

