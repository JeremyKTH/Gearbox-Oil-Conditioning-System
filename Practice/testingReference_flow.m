clc
clear all

%Initial Notes...
%   ~12 min to heat up 

%% Desired Values (INPUT HERE)

inputRPM = 1500;    %Speed of 
torqueInc = 500;    %Amount to increment torque each loop
maxTorque = 1500;

equilFlow = 1750;     %flow rate rpm (50% - 11 L/min and 1750 rpm)
maxFlow = 3500;    %3500
minFlow = 0;     %0
flowInc = 350;      %11

% Time variables:

heatTime = 20; %time to heat-up [min]
startTime = 20; % when to start the square wave[min]
coolDownTime = -15; % time to cooldown [min]
sampleRate = .1; % desired rate for sending reference [sec]

loadTime = 5;
equilTime = 15;
period = loadTime + equilTime;

altPeriod = 10; %time before switching back the other way


%% Conversions (min --> sec)

heatTime = heatTime*60/sampleRate; %sec
startTime = startTime*60/sampleRate;
coolDownTime = coolDownTime*60/sampleRate;
altPeriod = altPeriod*60/sampleRate;

loadTime = loadTime*60/sampleRate;
equilTime = equilTime*60/sampleRate;
period = period*60/sampleRate;

%% Setting up sim time

counterMax = (maxFlow-equilFlow)/flowInc;
%altStart = (startTime + period*(counterMax))*(60/sampleRate);

totalTime = (startTime+coolDownTime+altPeriod+(period*counterMax)*2)*sampleRate; %[sec]
halfTime = startTime + (period*counterMax) + (altPeriod);

time = 0:sampleRate:totalTime;

%% Generating reference 

for i=1:size(time,2)
    if i <= heatTime
        torque(i) = 0;
        flow(i) = equilFlow;
    elseif i >= (startTime + counterMax*period + 1)
        torque(i) = 0;
        flow(i) = equilFlow;
    else
        torque(i) = 0;
        flow(i) = equilFlow;
    end
end

for n = 0:counterMax-1
    for i = (startTime + period*n):(startTime + (period*(n+1))+1)
        if i < (startTime + period*(n+1) - equilTime)
            torque(i) = torqueInc;
            flow(i) = equilFlow + flowInc*(n+1);
        else
            torque(i) = 0;
        end
    end
    for i = halfTime:halfTime+altPeriod+(period*(n-1))
        torque(i) = torque(i-altPeriod-(period*counterMax-1));
        flow(i) = -flow(i-altPeriod-(period*counterMax-1)) + equilFlow*2;
    end
end
    
figure(1)
clf 
hold on

yyaxis right
plot(time./60, torque,'r', 'DisplayName', 'Torque Load [Nm]')
yticks([-1500:500:1500]);
ylim([-1750 1750]);
ylabel('Torque Load and Gearbox RPM')
yline(inputRPM , '--m', 'DisplayName', 'Gearbox RPM')

yyaxis left
plot(time./60, flow, 'DisplayName', 'Oil Pump RPM')
yticks([0:350:3500]);
ylim([-350 3850]);
ylabel('Oil Pump RPM')

%xticks([0:10:totalTime/60])
xlim([0 totalTime/60]);
%xlim([0 totalTime])
%yline(inputRPM , '--r', 'DisplayName', 'Gearbox RPM')
%yline(flowRPM , '--m', 'DisplayName', 'Oil Pump RPM (~11 L/min)')
title('Controller Test - FLow Disturbance')
xlabel('Time [min]')

legend()
%m = [time; torque]';
%writematrix(m, 'systemID_ReferenceInput.csv');
hold off