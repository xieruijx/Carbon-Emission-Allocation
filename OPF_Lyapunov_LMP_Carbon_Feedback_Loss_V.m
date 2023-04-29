clear
close all

tic

M = 50;
load Grid.mat
load('CCarbon.mat','CCarbon');
delta = 0.0001;
Cgen_k1 = Cgen_k+delta*CCarbon;
kappa_carbon = 0.05;
carbon_l = [0.00;0.00];
carbon_u = [0.31;1.45];
lambda_u = [0.55;1.78]+carbon_u*kappa_carbon;
lambda_l = [0.00;0.47]+carbon_l*kappa_carbon;
% Optimal
eta = (lambda_u*eta_c*eta_d*w_h-lambda_l*w_l).*ESe./(lambda_u*eta_c*eta_d-lambda_l);
V = (eta-w_l*ESe)./lambda_u/eta_d;

V(1) = 0.1*V(1);
eta(1) = (w_l*ESe(1)+V(1)*lambda_u(1)*eta_d+w_h*ESe(1)+V(1)*lambda_l(1)/eta_c)/2;

psSet = zeros(Nes,Num_T);
psfSet = zeros(Nes,Num_T);
pgSet = zeros(Ngen,Num_T);
prSet = zeros(NR,Num_T);
lambdaSet = zeros(Nbus,Num_T);
qsSet = zeros(Nes,Num_T);
CPriceSet = zeros(Nbus,Num_T);
% es0 = (w_l+w_h)*ESe/2;
es0 = w_h*ESe;
qs = es0-eta;
carbonES = carbon_l;
for t = 1:Num_T
    t
    PSmin = max(-ESp,(eta_c*qs+V.*lambda_l)/eta_c/eta_c/Delta_t);
    PSmin(qs>=-V.*lambda_l/eta_c) = 0;
    PSmax = min(ESp,(qs*eta_d+V.*lambda_u*eta_d*eta_d)/Delta_t);
    PSmax(qs<=-V.*lambda_u*eta_d) = 0;
    
    PSset = zeros(M,Nes);
    for i = 1:Nes
        PSset(:,i) = (linspace(PSmin(i),PSmax(i),M))';
    end
    fset = zeros(size(PSset));
    for i = 1:Nes
        fset(PSset(:,i)>0,i) = PSset(PSset(:,i)>0,i).*((PSset(PSset(:,i)>0,i)*Delta_t/2-qs(i)*eta_d)/V(i)/eta_d/eta_d-carbonES(i)*kappa_carbon);
        fset(PSset(:,i)<0,i) = PSset(PSset(:,i)<0,i).*((PSset(PSset(:,i)<0,i)*eta_c*eta_c*Delta_t/2-qs(i)*eta_c)/V(i)-carbonES(i)*kappa_carbon);
    end

    ps = sdpvar(Nes,1); % dischange-charge
    psf = sdpvar(Nes,1);
    pg = sdpvar(Ngen,1);
    pr = sdpvar(NR,1);

    obj_opf = sum(Cgen_k1'*pg*Delta_t)+sum(psf)*Delta_t;
    con1_opf = (1-losskappa)'*(Igen*pg+IR*pr+Ies*ps-PD0(:,t)) == 0;
    con2_opfa = PTDF_l*(Igen*pg+IR*pr+Ies*ps-PD0(:,t)) <= Sbranch_l;
    con2_opfb = -PTDF_l*(Igen*pg+IR*pr+Ies*ps-PD0(:,t)) <= Sbranch_l;
    con3_opf = [PMINgen <= pg <= PMAXgen,0 <= pr <= PR(:,t)];
    con4_opf = PSmin<=ps<=PSmax;
    for i = 1:Nes
        con4_opf = [con4_opf,psf(i)>=fset(1:M-1,i)+(fset(2:M,i)-fset(1:M-1,i))./(PSset(2:M,i)-PSset(1:M-1,i)).*(ps(i)-PSset(1:M-1,i))];
    end
    con_opf = [con1_opf,con2_opfa,con2_opfb,con3_opf,con4_opf];
    optimize(con_opf,obj_opf,sdpsettings('verbose',0));
    
    lambdaSet(:,t) = (-dual(con1_opf)*(1-losskappa)-PTDF_l'*dual(con2_opfa)+PTDF_l'*dual(con2_opfb))/Delta_t;
    psSet(:,t) = value(ps);
    psfSet(:,t) = value(psf);
    pgSet(:,t) = value(pg);
    prSet(:,t) = value(pr);
    qs = qs-eta_c*min(psSet(:,t),0)*Delta_t-max(psSet(:,t),0)*Delta_t/eta_d;
    qsSet(:,t) = qs;

    CPriceSet(:,t) = calCPrice_PTDF(PD0(:,t)-Ies*psSet(:,t),PR(:,t));
    carbonES = CPriceSet(Ines,t);
end

esSet = qsSet+eta*ones(1,Num_T);

figure;
plot(1:Num_T,lambdaSet(Ines(1),:),'LineWidth',2);
xlabel('Time (h)');
ylabel('LMP at Bus 15 ($/kWh)');
set(gca,'FontName','Times New Roman','FontSize',14);

figure;
plot(1:Num_T,lambdaSet(Ines(2),:),'LineWidth',2);
xlabel('Time (h)');
ylabel('LMP at Bus 18 ($/kWh)');
set(gca,'FontName','Times New Roman','FontSize',14);

figure;
plot(1:Num_T,CPriceSet(Ines(1),:),'LineWidth',2);
xlabel('Time (h)');
ylabel('Emission price at Bus 15 (kgCO_2/kWh)');
set(gca,'FontName','Times New Roman','FontSize',14);

figure;
plot(1:Num_T,CPriceSet(Ines(2),:),'LineWidth',2);
xlabel('Time (h)');
ylabel('Emission price at Bus 18 (kgCO_2/kWh)');
set(gca,'FontName','Times New Roman','FontSize',14);

figure;
plot(1:Nbus,CPriceSet(:,50),'LineWidth',2);
xlabel('Bus');
ylabel('Emission price (kgCO_2/kWh)');
set(gca,'FontName','Times New Roman','FontSize',14);

income = zeros(Num_T,2);
income(1,1) = (lambdaSet(Ines(1),1)+kappa_carbon*CPriceSet(Ines(1),1))*psSet(1,1);
for t = 2:Num_T
    income(t,1) = income(t-1,1)+(lambdaSet(Ines(1),t)+kappa_carbon*CPriceSet(Ines(1),t))*psSet(1,t);
end
income(1,2) = (lambdaSet(Ines(2),1)+kappa_carbon*CPriceSet(Ines(2),1))*psSet(2,1);
for t = 2:Num_T
    income(t,2) = income(t-1,2)+(lambdaSet(Ines(2),t)+kappa_carbon*CPriceSet(Ines(2),t))*psSet(2,t);
end
ab1 = polyfit((1:Num_T)',income(:,1),1);
ab2 = polyfit((1:Num_T)',income(:,2),1);

figure;
plot(1:Num_T,income(:,1),'LineWidth',2);
xlabel('Time (h)');
ylabel('Revenue of ES at bus 15 (10^4 $)');
set(gca,'FontName','Times New Roman','FontSize',14);

figure;
plot(1:Num_T,income(:,2),'LineWidth',2);
xlabel('Time (h)');
ylabel('Revenue of ES at bus 18 (10^4 $)');
set(gca,'FontName','Times New Roman','FontSize',14);

CCostSet_PD = CPriceSet.*PD0(:,1:Num_T);
CCostSet_ES = -psSet.*((Ies')*CPriceSet);
CCostAccum = zeros(Num_T+1,2);
for t = 1:Num_T
    CCostAccum(t+1,1) = CCostAccum(t,1)+CCostSet_ES(1,t);
end
for t = 1:Num_T
    CCostAccum(t+1,2) = CCostAccum(t,2)+CCostSet_ES(2,t);
end
ab1 = polyfit((1:Num_T+1)',CCostAccum(:,1),1);
ab2 = polyfit((1:Num_T+1)',CCostAccum(:,2),1);

figure;
plot(0:Num_T,CCostAccum(:,1),'LineWidth',2);
xlabel('Time (h)');
ylabel('Cumulative emission of ES at bus 15 (10^5 kgCO_2');
set(gca,'FontName','Times New Roman','FontSize',14);

figure;
plot(0:Num_T,CCostAccum(:,2),'LineWidth',2);
xlabel('Time (h)');
ylabel('Cumulative emission of ES at bus 18 (10^5 kgCO_2');
set(gca,'FontName','Times New Roman','FontSize',14);

CCostAccumSys = zeros(Num_T+1,1);
for t = 1:Num_T
    CCostAccumSys(t+1) = CCostAccumSys(t)+2*sum(CCostSet_PD(:,t))+2*sum(CCostSet_ES(:,t));
end
figure;
plot(0:Num_T,CCostAccumSys,'LineWidth',2);
xlabel('Time (h)');
ylabel('Cumulative system emission (10^5 kgCO_2');
set(gca,'FontName','Times New Roman','FontSize',14);

sum(CCarbon'*pgSet*Delta_t)/Num_T
toc
save('Result_Lyapunov_ES_Feedback_Loss.mat','PD0','psSet','pgSet','lambdaSet','PR','Num_T','Delta_t','kappa_carbon','CPriceSet','income','CCostAccum','CCostAccumSys');