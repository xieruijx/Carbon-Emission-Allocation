clear
close all
load Grid.mat

tNum = 50;
tNum = 13;

%% Proposed

[CPrice,~] = calCPrice_PTDF_num(PD0(:,tNum),PR(:,tNum));

figure;
plot(1:Nbus,CPrice,'LineWidth',2);
xlabel('Bus');
ylabel('Emission price (kgCO_2/kWh)');
set(gca,'FontName','Times New Roman','FontSize',14);

%% CEF

load('CCarbon.mat','CCarbon');
delta = 0.0001;
Cgen_k1 = Cgen_k+delta*CCarbon;

pg = sdpvar(Ngen,1);
pr = sdpvar(NR,1);
pl = sdpvar(Nbranch,1);
theta = sdpvar(Nbus,1);

obj = Cgen_k1'*pg;
con1 = pg >= PMINgen;
con2 = pg <= PMAXgen;
con3 = pl >= -Sbranch;
con4 = pl <= Sbranch;
con5 = ITo'*pl+Igen*pg+IR*pr-IFrom'*pl == PD0(:,tNum);
con6 = BR_X.*pl == (IFrom-ITo)*theta;
con7 = [theta(1) == 0, 0 <= pr <= PR(:,tNum)];
con = [con1,con2,con3,con4,con5,con6,con7];

optimize(con,obj);

Pl = value(pl);
Pg = [value(pg);value(pr)];
CCarbon = [CCarbon;0;0];
Igen = [Igen,IR];

for l = 1:Nbranch
    if Pl(l) < 0
        temp = Ibranch(l,1);
        Ibranch(l,1) = Ibranch(l,2);
        Ibranch(l,2) = temp;
        ITo(l,:) = zeros(1,Nbus);
        ITo(l,Ibranch(l,2)) = 1;
        IFrom(l,:) = zeros(1,Nbus);
        IFrom(l,Ibranch(l,1)) = 1;
        Pl(l) = -Pl(l);
    end
end
assert(all(Pl>=0))
PlMatrix = zeros(Nbus,Nbus);
for l = 1:Nbranch
    PlMatrix(Ibranch(l,2),Ibranch(l,1)) = PlMatrix(Ibranch(l,2),Ibranch(l,1))+Pl(l);
end
ZIndex = Igen*Pg+ITo'*Pl == 0;
if sum(ZIndex) > 0
    Igenn = Igen(~ZIndex,:);
    ITon = ITo(:,~ZIndex);
    PlMatrix = PlMatrix(~ZIndex,~ZIndex);
    NZrho = (diag(Igenn*Pg+ITon'*Pl)-PlMatrix)\(Igenn*(CCarbon.*Pg));
else
    NZrho = (diag(Igen*Pg+ITo'*Pl)-PlMatrix)\(Igen*(CCarbon.*Pg));
end
rhoSet = zeros(Nbus,1);
rhoSet(~ZIndex) = NZrho;

figure;
plot(1:Nbus,rhoSet,'LineWidth',2);
xlabel('Bus');
ylabel('Emission price (kgCO_2/kWh)');
set(gca,'FontName','Times New Roman','FontSize',14);

figure;
bar([CPrice(PD0(:,tNum)>0),rhoSet(PD0(:,tNum)>0)/2]);
legend('Proposed','CEF');
xlabel('Bus');
ylabel('Emission price (kgCO_2/kWh)');
xticks(1:20);
xticklabels({'2','3','4','7','8','10','12','14','15','16','17','18','19','20','21','23','24','26','29','30'});
set(gca,'FontName','Times New Roman','FontSize',14);