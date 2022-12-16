%% Home of the Simulink Model
% 11/9/22
% main3 has some food info, but the paper does too

%%
clear; clc; format compact; close all

%% Variables for Simulink

n_drink = 1; % number of standard drinks taken 
W = 90; % Mass of the subject (kg)
t_end_min = 800; % Time of the simulation in minutes (min)

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
k_p = 0.39512;
k_i = 0.58368;
k_d = 0;

wC_out = sim("Proj_v3_wF_wC.slx",t_end);
t_cont = wC_out.simout.Time/60; % Convert s to mins
BEC_cont = wC_out.simout.Data/10; % Convert L to dL
tf_cont = wC_out.food.Time/60;
food_cont = wC_out.food.Data;
error = wC_out.error.Data;

woC_out = sim("Proj_v3_wF.slx",t_end);
t_no = woC_out.simout.Time/60; % Convert s to mins
BEC_no = woC_out.simout.Data/10; % Convert L to dL
tf_no = woC_out.food.Time/60;
food_no = woC_out.food.Data;
M_GV_no = woC_out.M_GV.Data;

figure;
subplot(311);plot(t_cont,BEC_cont,"o",t_no,BEC_no,"o");xlabel("Time (min)");ylabel("BEC (mg/dL)");title("Comparing Controller vs. No Controller");legend("Controller","N/A");xlim([0 4])
subplot(312);plot(tf_cont,food_cont,"o",tf_no,food_no,"o");xlabel("Time (min)");ylabel("Food");legend("Controller","N/A");xlim([0 4])
subplot(313);plot(t_cont,error,"o");title("Error");xlabel("Time (min)");xlim([0 4])
close

figure;
subplot(311);plot(t_cont,BEC_cont,t_no,BEC_no);xlabel("Time (min)");ylabel("BEC (mg/dL)");title("Comparing Controller vs. No Controller");legend("Controller","N/A");
subplot(312);plot(tf_cont,food_cont,tf_no,food_no);xlabel("Time (min)");ylabel("Food");legend("Controller","N/A");
subplot(313);plot(t_cont,error);title("Error");xlabel("Time (min)");

% figure;
% s = tf([.1736],[1 0])
% rlocus(s)
%}

%% tfest()
% https://www.mathworks.com/help/ident/ref/iddata.html
% https://www.mathworks.com/help/ident/ug/data-types-in-system-identification-toolbox.html#DataTypesInSystemIdentificationToolboxExample-8
% https://www.mathworks.com/help/ident/ref/tfest.html#responsive_offcanvas
% https://www.mathworks.com/help/ident/ref/tfestoptions.html

woC_out = sim("Proj_v3_wF.slx",t_end);
t_no = woC_out.simout.Time/60; % Convert s to mins
BEC_no = woC_out.simout.Data/10; % Convert L to dL
M_GV_no = woC_out.M_GV.Data;

measTime = t_no;
% input = M_GV_no;
% output = BEC_no;
input = ones([1 numel(measTime)])*r;
input = input';
output = M_GV_no;

% load y_min_yHat.mat
% measTime = t_no(1:numel(y_min_yHat));
% input = M_GV_no(1:numel(y_min_yHat));
% output = y_min_yHat;



% tt = iddata(output,input,"InputName",["M_GV"],"InputUnit",["mg/s"],"OutputName",["BEC"],"OutputUnit",["mg/dL"],"SamplingInstants",measTime,"TimeUnit",["minutes"]);
tt = iddata(output,input,"InputName",["r"],"InputUnit",["mg/s"],"OutputName",["M_GV"],"OutputUnit",["mg/dL"],"SamplingInstants",measTime,"TimeUnit",["minutes"]);

%{
row_n = 3;
for nz = 1:5
    for np = 1:5
        try
            tf_est = tfest(tt,np,nz,"InputName",["M_GV"],"OutputName",["BEC"]);
            WriteToExcel("Test.xlsx",tf_est,np,nz,row_n)
            row_n = row_n + 1;
        catch
            disp("Ignored!")
        end
    end
end
%}

% est_tf = tfest([measTime input1 input2 output],np);
% estimate_tf = tfest([t_no BEC_no],np)
% Look for positive % closest to 100%

function WriteToExcel(f_name,tf,np,nz,row_n)
row_n = num2str(row_n);

writematrix(tf.Report.Fit.FitPercent,f_name,"Range",append("A",row_n))
writematrix(nz,f_name,"Range",append("B",row_n))
writematrix(np,f_name,"Range",append("C",row_n))

end
