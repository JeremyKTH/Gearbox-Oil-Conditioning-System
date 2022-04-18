% Open Loop --------------------
% load 'ARMAX_open_sine_Err.mat'
% load 'ARMAX_open_square_Err.mat'
% load 'TF_open_sine_Err.mat'
% load 'TF_open_square_Err.mat'
% load 'SS_open_sine_Err.mat'
load 'SS_open_square_Err.mat'

% Closed loop ------------------
% load 'ARMAX_sine_Err.mat'
% load 'ARMAX_square_Err.mat'
% load 'TF_sine_Err.mat'
% load 'TF_square_Err.mat'
% load 'SS_sine_Err.mat'
% load 'SS_square_Err.mat'

sum = 0;

for i = 1: length(data.Data)
    x = data.Data(i);
    sum = sum + x^2;
end

NRMSE = sqrt(sum/length(data.Data))