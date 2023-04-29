clear
clc
close all

% mpc = loadcase('case30pwl.m');
mpc = loadcase('case118.m');
PD = mpc.bus(:,3) / mpc.baseMVA;
day = 28;
Delta_t = 1;
sigma = 0.1;
Num_T = day/Delta_t*24;
daily = ([0.56,0.53,0.54,0.51,0.50,0.50,...
    0.53,0.59,0.62,0.64,0.64,0.63,...
    0.65,0.62,0.61,0.61,0.59,0.57,...
    0.59,0.65,0.63,0.58,0.55,0.59]-0.4)/0.25;
base = interp1(0:24,[daily,daily(1)],0:Delta_t:23.99);
PDdata = zeros(length(PD),Num_T);
for i = 1:day
    temp = PD*base;
    temp = temp.*(1+sigma*randn(length(PD),24/Delta_t));
    PDdata(:,1+(i-1)*24/Delta_t:i*24/Delta_t) = temp;
end
PDdata = PDdata*1.0;
load('PVdata.mat','PV');
PVdata = interp1(0:24*day,PV(1:24*day+1),0:Delta_t:24*day-0.01);
load('Winddata.mat','Wind');
Winddata = interp1(0:24*day,Wind(1:24*day+1),0:Delta_t:24*day-0.01);
PDu = [3.70;3.70;3.74;3.74;3.73;3.69;3.73;3.72;3.80;3.84;
    3.80;3.82;3.78;3.83;3.81;3.83;3.84;3.84;3.84;3.84;
    3.86;3.83;3.86;3.98;4.26;4.26;3.30;3.69;3.34;3.34];
% save('Data30.mat','PDdata','PVdata','Winddata','PDu');
save('Data118.mat','PDdata','PVdata','Winddata','PDu');