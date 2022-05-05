clc
clear all

%Initial Notes...
%   ~12 min to heat up 

%% Desired Values (INPUT HERE)

inputRPM = 1500;    %Speed of 
torqueInc = 250;    %Amount to increment torque each loop
flowRPM = 1750;     %flow rate rpm (50% - 11 L/min and 1750 rpm)
maxTorque = 1250;

% Time variables:
%equil = 55; % equilibrium point for valve openeing [%]
heatTime = 10; %time to heat-up [min]
coolDownTime = 5; % time to cooldown [min]
sampleRate = .1; % desired rate for sending reference [sec]
startTime = 10; % when to start the square wave[min]
period = 20; % period of increment [min]

loadTime = 5;
equilTime = 15;
period = loadTime + equilTime;

%% Conversions (min --> sec)

heatTime = heatTime*60/sampleRate; %sec
startTime = startTime*60/sampleRate;
coolDownTime = coolDownTime*60/sampleRate;

loadTime = loadTime*60/sampleRate;
equilTime = equilTime*60/sampleRate;
period = period*60/sampleRate;

%% Setting up sim time

counterMax = maxTorque/torqueInc;

totalTime = (startTime+coolDownTime+(period*counterMax))*sampleRate; %[sec]

time = 0:sampleRate:totalTime;

%% Generating reference 

for i=1:size(time,2)
    if i <= heatTime
        torque(i) = 0;
    elseif i >= (startTime + counterMax*period + 1)
        torque(i) = 0;
    else
        torque(i) = 0;
    end
end

for n = 0:counterMax-1
    for i = (startTime + period*n):(startTime + (period*(n+1))+1)
        if i < (startTime + period*(n+1) - equilTime)
            torque(i) = torqueInc*(n+1);
        else
            torque(i) = 0;
        end
    end
end

figure(1)
clf 
plot(time./60, torque, 'DisplayName', 'Torque Load [Nm]')
%plot(time, inputRPM )
yticks([-250:250:2250]);
ylim([-250 2250]);
xticks([0:10:totalTime/60])
xlim([0 totalTime/60]);
%xlim([0 totalTime])
yline(inputRPM , '--r', 'DisplayName', 'Gearbox RPM')
yline(flowRPM , '--m', 'DisplayName', 'Oil Pump RPM (~11 L/min)')
title('Controller Test - Torque Disturbance')
xlabel('Time [min]')
ylabel('Torque and RPM')
legend()
%m = [time; torque]';
%writematrix(m, 'systemID_ReferenceInput.csv');
