%% Best +/-40% TF Model
b = 0.001273;
a = -0.9965;
B = [b];         % B
A = [1, a];  % A
Ts = 0.4;

Gp = tf(B, A, Ts);

poles = pole(Gp);  % 0.9965
zeros = zero(Gp);  % No zeros
%pzmap(Gp)
%step(Gp)
%% Choose Poles (w_m) (w_o)
%--- A_m ------
w_m = 1;
%--- A_o ------
w_o = 1;

p1 = exp(-w_m*Ts);
p0 = exp(-w_o*Ts);

fprintf('The original disc. poles are: 0.9965 \n');
fprintf('The chosen cont. poles are: %d and %d \n', -p1, -p0);
chosen_disc_poles = ['The chosen disc. poles are: ', num2str(p1), ' and ', num2str(p0)];
disp(chosen_disc_poles)

%% Diophantine Eqn A_cl = AR + BS
syms P I z
A_cl = (z+a)*((z-1)*P + I*Ts);
A_d = (z-p1)*(z-p0);
A_cl_c = fliplr(coeffs(A_cl, z)); % retreive coeficients
A_d_c = fliplr(coeffs(A_d, z));   % retreive coeficients

equ1 = A_cl_c(2) == A_d_c(2); % z^1
equ2 = A_cl_c(3) == A_d_c(3); % z^0

sol = solve([equ1, equ2], [P, I]);

P = double(sol.P)   % 0.8929
I = double(sol.I)   % 1.1049

%% Gc
z = tf('z', Ts);
Gc = P + I*Ts/(z-1);

%% Gff and Gyr
t_o = (1-p1)/(b);
A_o = [1, -p0];
T = t_o*A_o;
R = [1, -1];
Gff = tf(T, R, Ts);

%% Entire closed loop check
Gyr = Gff*Gp/(1+Gc*Gp)
pole(Gyr)
zero(Gyr)

figure(1)
pzmap(Gyr)

figure(2)
step(Gyr)