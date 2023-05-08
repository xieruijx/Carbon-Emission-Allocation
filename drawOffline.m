clear
close all

load Grid.mat
load Result_Lyapunov_ES_Feedback_Loss_V1_k0.05.mat

Num_T = 672;
income = income(1:Num_T,:);

Ines = 15;
% ESe(1) = ESe(1)*2;

gammaSet = lambdaSet(Ines,:)+kappa_carbon*CPriceSet(Ines,:);
gammaSet = gammaSet(:,1:Num_T);
es0 = w_h*ESe(1);

%% Offline

psc = sdpvar(Num_T,1);
psd = sdpvar(Num_T,1);
es = sdpvar(Num_T,1);
obj = gammaSet*(psc-psd);
con_ES = [0<=psc<=ESp(1),0<=psd<=ESp(1),es==[es0;es(1:end-1)]+psc*Delta_t*eta_c-psd*Delta_t/eta_d,w_l*ESe(1)<=es<=w_h*ESe(1)];
optimize(con_ES,obj);

PSC = value(psc);
PSD = value(psd);
PS = PSD-PSC;

Revenue_Off_Theo = zeros(Num_T,1);
for t = 1:Num_T
    Revenue_Off_Theo(t) = gammaSet(1:t)*PS(1:t);
end

%% Traditional drift-plus-penalty

load Grid.mat
load Result_Lyapunov_ES_Feedback_Loss_V1_k0.05.mat

Num_T = 672;
income = income(1:Num_T,:);

Ines = 15;

gammaSet = lambdaSet(Ines,:)+kappa_carbon*CPriceSet(Ines,:);
gammaSet = gammaSet(:,1:Num_T);
es0 = w_h*ESe(1);

carbon_l = [0.00;0.00];
carbon_u = [0.31;1.45];
lambda_u = [0.55;1.78]+carbon_u*kappa_carbon;
lambda_l = [0.00;0.45]+carbon_l*kappa_carbon;

V = (w_h*ESe(1)-w_l*ESe(1)-ESp(1)*eta_c*Delta_t-ESp(1)*Delta_t/eta_d)/(lambda_u(1)*eta_d-lambda_l(1)/eta_c);
eta = w_h*ESe(1)+V*lambda_l(1)/eta_c-ESp(1)*eta_c*Delta_t;
qs = w_h*ESe(1)-eta;

pscSet = zeros(Num_T,1);
psdSet = zeros(Num_T,1);
qsSet = zeros(Num_T,1);
% Online
for t = 1:Num_T
    psc = 0;
    psd = 0;
    if qs <= -V*gammaSet(t)/eta_c
        psc = ESp(1);
    elseif qs >= -V*gammaSet(t)*eta_d
        psd = ESp(1);
    end
    
    pscSet(t) = psc;
    psdSet(t) = psd;
    qs = value(qs+eta_c*psc*Delta_t-psd*Delta_t/eta_d);
    qsSet(t) = qs;
end

Revenue_Tra_Theo = zeros(Num_T,1);
for t = 1:Num_T
    Revenue_Tra_Theo(t) = gammaSet(1:t)*(psdSet(1:t)-pscSet(1:t));
end

%% Traditional price bound

load Grid.mat
load Result_Lyapunov_ES_Feedback_Loss_V1_k0.05.mat

Ines = 15;
Num_T = 672;

gammaSet = lambdaSet(Ines,:)+kappa_carbon*CPriceSet(Ines,:);
gammaSet = gammaSet(:,1:Num_T);
es0 = w_h*ESe(1);

gamma_u = 0.5;
gamma_l = 0.3;

pscSet = zeros(Num_T,1);
psdSet = zeros(Num_T,1);
es = es0;
for t = 1: Num_T
    if gammaSet(t)>gamma_u
        psdSet(t) = min(ESp(1),(es-w_l*ESe(1))/Delta_t*eta_d);
    elseif gammaSet(t)<gamma_l
        pscSet(t) = min(ESp(1),(w_h*ESe(1)-es)/Delta_t/eta_c);
    end
    es = es+pscSet(t)*eta_c*Delta_t-psdSet(t)*Delta_t/eta_d;
end

Revenue_Bound_Theo = zeros(Num_T,1);
for t = 1:Num_T
    Revenue_Bound_Theo(t) = gammaSet(1:t)*(psdSet(1:t)-pscSet(1:t));
end

figure;
plot(0:Num_T,[0;income(:,1)],0:Num_T,[0;Revenue_Off_Theo],0:Num_T,[0;Revenue_Tra_Theo],0:Num_T,[0;Revenue_Bound_Theo],'LineWidth',2);
legend('Proposed','Offline','Traditional real-time','Simple','location','northwest');
xlabel('Time (h)');
ylabel('Revenue of ES 15 (10^4 $)');
set(gca,'FontName','Times New Roman','FontSize',14);

figure;
plot(gammaSet);