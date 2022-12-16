%% Home of the Simulink Model
% 11/9/22
% main3 has some food info, but the paper does too

%%
clear; clc; format compact; close all

%% Variables for Simulink

n_drink = 1; % number of standard drinks taken 
W = 90; % Mass of the subject (kg)
t_end_min = 1000; % Time of the simulation in minutes (min)

t_end = t_end_min*60; % End time of the simulation (s)
drink_m = 14; % mass of EtOH per standard drink (g)
dose = drink_m*1000*n_drink; % mass of ethanol (mg) in drink
L = 60; % Length of the drink (s)
S = 0; % Start time of the drink (s)
r = dose/L;
b = 0.0007363716; % A constant
s = 1.5; % A constant for food
m = 0.75; % A constant for food
delta = 3.1e-4; % Digestive rate (1/s)
E = S + L; % The end time of the drink

rho_liver = 1.13; % g/mL
% Volumes
V_v = 0.064*W; % L
V_L = (((W*0.0148)*1000)/rho_liver)/1000; % L
V_T = (W/1.01)*(1-0.064-0.12-0.0148); % L

% V_L = 1.2863; % L
% V_v = 0.64; % L
% V_T = 13.404; % L

% Maximum motabolism rate
V_max = 3.0265; % mg/(L*min)

% Ks
k = 0.7931544;
K_m = 16.55418;

% Volume flow rate of aorta
R_A = 0.0705;

% Fraction of cardiac outputs
F_P = 0.75;
F_L = 0.26;

% PI Controller coefficients
k_p = 0.053032;
k_i = 0.00054959;
% k_p = 0.053032;
% k_i = 0.00054959;

%% Run the simulation

% wF_out = sim("Proj_v3_wF.slx",t_end); 
% wF_time = wF_out.simout.Time/60; % Converting seconds to minutes
% wF_BEC = wF_out.simout.Data/10; % Converting L to dL
% 
% out = sim("Proj_v3.slx",t_end);
% time = out.simout.Time/60; % Converting seconds to hours
% BEC = out.simout.Data/10; % Converting L to dL
%
% This will plot both models on top of each other
% figure;
% plot(wF_time,wF_BEC,time,BEC);legend("w/ F","w/out F");xlabel("Time (min)");ylabel("BEC (mg/dL)");grid on
%
% [max_out,i_max] = max(out.simout.Data); % Find the maximum value and the index
% out_time_max = out.simout.Time(i_max)/60; % Time of the peak
% out_BEC_max = out.simout.Data(i_max)/10; % Peak BEC
% fprintf("Peak of %4.4f mg/dL at %4.4f min\n",out_BEC_max,out_time_max)
% figure;plot(time,BEC);xlabel("Time (min)");ylabel("BEC (mg/dL)");title("One Drink");
%
% out = sim("Proj_v3_step.slx",t_end);
% time = out.simout.Time/60; % Converting seconds to mins
% BEC = out.simout.Data/10; % Converting L to dL
%
% arb_info = lsiminfo(BEC,time)
% step_info = stepinfo(BEC,time)
% figure;plot(time,BEC);xlabel("Time (min)");ylabel("BEC (mg/dL)");title("One Drink");
%
% figure;
% plot(out_time,out_BEC,out_time_max,out_BEC_max,"o");xlabel("Time (min)");ylabel("BEC (mg/dL)");title("Output");grid on


%% Tuning the Controller
%{
wC_out = sim("Proj_v3_wF_wC.slx",t_end);
t_cont = wC_out.simout.Time/60; % Convert s to mins
BEC_cont = wC_out.simout.Data/10; % Convert L to dL
tf_cont = wC_out.food.Time/60;
food_cont = wC_out.food.Data;
error = wC_out.error.Data;
[max_F,i_F] = max(food_cont);
[max_w,i_w_max] = max(BEC_cont);
fprintf("Maximum food is %4.2f grams at %4.2f minutes.\n",max_F,tf_cont(i_F))
fprintf("Maximum BEC (with PI) is %4.2f mg/dL at %4.2f minutes.\n",max_w,tf_cont(i_w_max))

woC_out = sim("Proj_v3_wF.slx",t_end);
t_no = woC_out.simout.Time/60; % Convert s to mins
BEC_no = woC_out.simout.Data/10; % Convert L to dL
tf_no = woC_out.food.Time/60;
food_no = woC_out.food.Data;
M_GV_no = woC_out.M_GV.Data;
[max_wo,i_wo_max] = max(BEC_no);
fprintf("Maximum BEC (without PI) is %4.2f mg/dL at %4.2f minutes.\n",max_wo,tf_cont(i_wo_max))

figure;
subplot(311);plot(t_cont,food_cont);xlabel("Time (min)");ylabel("Food (g)")

% figure;
% subplot(311);plot(t_cont,BEC_cont,"o",t_no,BEC_no,"o");xlabel("Time (min)");ylabel("BEC (mg/dL)");title("Comparing Controller vs. No Controller");legend("Controller","N/A");xlim([0 4])
% subplot(312);plot(tf_cont,food_cont,"o",tf_no,food_no,"o");xlabel("Time (min)");ylabel("Food");legend("Controller","N/A");xlim([0 4])
% subplot(313);plot(t_cont,error,"o");title("Error");xlabel("Time (min)");xlim([0 4])
% figure;
% subplot(311);plot(t_cont,BEC_cont,t_no,BEC_no);xlabel("Time (min)");ylabel("BEC (mg/dL)");title("Comparing Controller vs. No Controller");legend("Controller","N/A");
% subplot(312);plot(tf_cont,food_cont,tf_no,food_no);xlabel("Time (min)");ylabel("Food");legend("Controller","N/A");
% subplot(313);plot(t_cont,error);title("Error");xlabel("Time (min)");

% figure;
% plot(t_no,BEC_no,t_cont,BEC_cont);xlabel("Time (min)");ylabel("BEC (mg/dL)");legend("Plant","PI Compensated");

%}

%% Brute force
%{
% https://www.mathworks.com/help/ident/ug/model-quality-metrics.html
% https://www.mathworks.com/help/ident/ref/tfest.html
% load Static_Input.mat

% num = [0.2455 0.00663];
% den = [1 1.341 0.3551 0.002696];

%{
% zeros = [-0.027];
% poles = [-0.9823 -0.3509 -0.0078];
% gains = [1];
zeros_f = [-0.027];
% poles = [-12 -1 -0.00018];
poles_f = [-12 -6 -0.02];
gains_f = [1];
[num_f,den_f] = zp2tf(zeros_f,poles_f,gains_f);

zeros_s = [];
% poles = [-0.0004];
poles_s = [-0.0004];
gains_s = [1];
[num_s,den_s] = zp2tf(zeros_s,poles_s,gains_s);

tf_f = poly2sym(num_f)/poly2sym(den_f);
tf_s = poly2sym(num_s)/poly2sym(den_s);
[num,den] = numden(simplify(tf_f + tf_s));
num = sym2poly(num);
den = sym2poly(den);
%}
num = [0.004784 -1.903e-10];
den = [1 2.662 0.2294 0.005025];

out = sim("Copy_of_Proj_v3_wF_TF.slx",t_end);
t = out.simout.Time/60; % Convert s to mins
BEC_big = out.simout.Data/10; % Convert L to dL
BEC_tf = out.tf_out.Data;

figure;
plot(t,BEC_big,t,BEC_tf);xlabel("Time (min)");ylabel("BEC (mg/dL)");legend("Model","TF");title("Comparing Plant and System Transfer Function")

% scaler = (max(BEC_big)/max(BEC_tf));
% BEC_tf = BEC_tf*scaler;

% figure;
% plot(t,BEC_big,t,BEC_tf);xlabel("Time (min)");ylabel("BEC (mg/dL)");legend("Model","TF");title("Comparing Plant and System Transfer Function")


y_min_yHat = BEC_big - BEC_tf;
% clearvars -except y_min_yHat
% save("y_min_yHat")

% RMSE is between 0 and 1, and 0 is perfect
RMSE = sqrt(sum((BEC_big-BEC_tf).^2)/numel(BEC_tf))
%}

%% Plotting the arguments of H(x)

%{
out = sim("Proj_v3_wF_Copy.slx",t_end);
t = out.simout.Time/60; % Convert s to mins
BEC = out.simout.Data/10; % Convert L to dL
input = out.M_GV.Data;

HL = out.M_HL_arg.Data;
LH = out.M_LH_arg.Data;
PL = out.M_PL_arg.Data;
TV = out.M_TV_arg.Data;
VT = out.M_VT_arg.Data;
alph = out.alph.Data;
uL_vL = out.uL_vL.Data;

figure;
subplot(311);plot(t,HL,t,LH,t,PL,t,TV,t,VT);legend("HL","LH","PL","TV","VT");ylim([0 50])
subplot(312);plot(t,HL,t,LH,t,PL,t,TV,t,VT);legend("HL","LH","PL","TV","VT");ylim([-25 0])
subplot(313);plot(t,alph,t,uL_vL);legend("alpha","mu_L/V_L")
% subplot(414);plot(t,BEC,t,input);legend("BEC","M_G_V");
% close
%}

%% 2 TF estimate models

%{
transient = sim("Proj_v3_wF_Transient.slx",t_end);
t_tran = transient.simout.Time/60; % Convert s to mins
BEC_tran = transient.simout.Data/10; % Convert L to dL
tf_tran = transient.tf_out.Data;

steady = sim("Proj_v3_wF_Steady.slx",t_end);
t_stea = steady.simout.Time/60; % Convert s to mins
BEC_stea = steady.simout.Data/10; % Convert L to dL
tf_stea = steady.tf_out.Data;

figure;
subplot(211);plot(t_tran,BEC_tran,t_tran,tf_tran);legend("Model","TF");title("Short-Term Response");xlabel("Time (min)");ylabel("BEC (mg/dL)")
subplot(212);plot(t_stea,BEC_stea,t_stea,tf_stea);legend("Model","TF");title("Long-Term Response");xlabel("Time (min)");ylabel("BEC (mg/dL)")
%}


%% Running Standard Responses

num = [1770 33610 127910 2540];
den = [2.5e6 4.5051e7 180918020 3672360 1440];

data = sim("Proj_v3_step.slx",t_end);

t = data.simout.Time/60; % s -> min
real_BEC = data.simout.Data/10; % L -> dL
tf_BEC = data.tf_out.Data;

figure;
plot(t,real_BEC,t,tf_BEC);xlabel("Time (min)");ylabel("BEC (mg/dL)");legend("Model","TF Est.")


