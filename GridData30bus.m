%% mpc data
clear
clc
close all

mpc = loadcase('case30pwl.m');
[Ngen,~] = size(mpc.gen); % number of generators
Cgen_k = mpc.gencost(:,12)./mpc.gencost(:,11)*mpc.baseMVA/1e4; % 1e4 $/p.u.
Ingen = mpc.gen(:,1);
Nbus = 30; % number of buses
Iref = 1; % reference bus
PD = mpc.bus(:,3) / mpc.baseMVA; % demand
Igen = zeros(Nbus, Ngen); % bus of generator
for i = 1: Ngen
    Igen(mpc.gen(i, 1), i) = 1;
end
PMAXgen = mpc.gen(1: Ngen, 9) / mpc.baseMVA;
PMINgen = mpc.gen(1: Ngen, 10) / mpc.baseMVA;

Ibranch = mpc.branch(:, 1: 2); % branch: from bus, to bus
[Nbranch, ~] = size(Ibranch);
BR_R = mpc.branch(:,3);
BR_X = mpc.branch(:,4);
BR_G = real(1./(BR_R+1i*BR_X));
BR_B = imag(1./(BR_R+1i*BR_X));
Sbranch = mpc.branch(:,6)/mpc.baseMVA; % branch capacity
IFrom = zeros(Nbranch, Nbus);
ITo = zeros(Nbranch, Nbus);
for i = 1: Nbranch
    IFrom(i, Ibranch(i, 1)) = 1;
    ITo(i, Ibranch(i, 2)) = 1;
end
PTDF = zeros(Nbranch+Nbus-1,Nbus);
Arel = IFrom(:,2:Nbus)-ITo(:,2:Nbus);
Mrel = [Arel',zeros(Nbus-1,Nbus-1);diag(BR_X),Arel];
for i = 2:Nbus
    erel = zeros(Nbranch+Nbus-1,1);
    erel(i-1) = 1;
    PTDF(:,i) = Mrel\erel;
end
PTDF = PTDF(1:Nbranch,:);
PTDF_l = [PTDF(1:12,:);PTDF(14:33,:);PTDF(35:41,:)];
Sbranch_l = [Sbranch(1:12);Sbranch(14:33);Sbranch(35:41)];


% Num_T = 96;
Num_T = 672;
Delta_t = 1;
eta_c = 0.95;
eta_d = 0.95;
w_l = 0.1;
w_h = 0.9;
Nes = 2; % number of energy storage
Ines = [15;18];
Ies = zeros(Nbus, Nes);
for i = 1:Nes
    Ies(Ines(i),i) = 1;
end
ESe = [0.4;0.2];
ESp = [0.04;0.04];

load('Data30.mat','PDdata','PVdata','Winddata','PDu');
PD0 = PDdata;
NR = 2;
InR = [6;15];
IR = zeros(Nbus,NR);
for i = 1:NR
    IR(InR(i),i) = 1;
end
PRMAX = [1;1];
PR = zeros(NR,Num_T);
PR(1,:) = PRMAX(1)*PVdata(1:Num_T);
PR(2,:) = PRMAX(2)*Winddata(1:Num_T);

losskappa = 0.005;
losskappa = -losskappa*ones(Nbus,1);
losskappa(Ingen) = -losskappa(Ingen);
losskappa(InR) = -losskappa(InR);


save Grid.mat

CCarbon = [0.9;0.8;0.8;0.2;0.3;0.3]; % kgCO_2/kWh
CCarbon = CCarbon(1:Ngen);


save('CCarbon.mat','CCarbon');