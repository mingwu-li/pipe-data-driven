% In this demo, we show how to modify the obtained ROM to impose control
% input.



close all;
clear all;
% fit SSM-based ROM without control in post-flutter region
transit_orb;
%% fit control matrix B in reduced dynamics
% Generate trajectories with control input. We do this by generating a
% smooth sequence of inputs, namely, through periodic actuation at various
% amplitudes and frequencies. The control is applied via base excitation at
% the clamped end.
amp  = [-0.3:0.05:-0.1 0.1:0.05:0.3];
omg  = 2:4:22;
namp = numel(amp);
nomg = numel(omg);
useq = cell(namp*nomg,1);
IC   = zeros(2*n,1);
xcData = cell(namp*nomg,2);
UData  = cell(namp*nomg,1);
numberPeriods = 100; numberPointsPerPeriod = 100;
nSamp = numberPeriods*numberPointsPerPeriod+1;
figure; hold on
for i=1:namp
    for j=1:nomg
        idx  = (i-1)*nomg+j;
        famp = amp(i); fom = omg(j);
        FForced = @(t,x,w) F(t,x) + [zeros(n,1); famp*(M\fExt)*cos(w*t)];
        opts = odeset('RelTol', 1e-4); 
        [t,x] = ode15s(@(t,x) FForced(t,x,fom),...
            0:dt:10, IC,opts);
        xcData{idx,1} = t.';
        xcData{idx,2} = (x(:,1:n)*phimid).';
        UData{idx,1}  = t.';
        UData{idx,2}  = famp*cos(fom*t.');
        plot(t,x(:,1:n)*phimid,'Linewidth',1);       
    end
end

%% Setup trajectories for fit
SSMDim = 2;
overEmbed = 0;
Yu = coordinatesEmbedding(xcData, SSMDim, 'OverEmbedding', overEmbed);
Udelay = coordinatesEmbedding(UData, SSMDim, 'OverEmbedding', overEmbed);

%% yuxData = cell(size(Yu,1),4);
for iTraj = 1:size(Yu, 1)
   yuxData{iTraj,1} = Yu{iTraj,1};
   yuxData{iTraj,2} = Yu{iTraj,2};
   yuxData{iTraj,3} = Udelay{iTraj,2}(1,:);
   yuxData{iTraj,4} = IMInfo.chart.map(Yu{iTraj,2});
end
%% Fit control matrices
[Br,regErrors] = fitControlMatrices(yuxData,RDInfo,10);
regErrors
%% Update models
Rauton = RDInfo.reducedDynamics.map;
R = @(t,x,u) Rauton(x) + Br*u;

%% comparision between updated model and direct solution
famp = 0.015; fom = 3;
FForced = @(t,x,w) F(t,x) + [zeros(n,1); famp*(M\fExt)*cos(w*t)];
opts = odeset('RelTol', 1e-4); 
[t,x] = ode15s(@(t,x) FForced(t,x,fom),...
    linspace(0,10,1000), IC,opts);
xmid = (x(:,1:n)*phimid).';

% Check dynamics
uFun = @(t) famp*cos(fom*t);
% x0 = IMInfo.chart.map(IC);
x0 = [0;0];
[tiPred,xiPred] = ode45(@(t,x) R(t,x,uFun(t)),t,x0); 
tiPred = transpose(tiPred); xiPred = transpose(xiPred); 
% figure; plot(tiPred,xiPred)


paramFun = IMInfo.parametrization.map;
zpred = paramFun(xiPred);

figure; plot(t,xmid,'r-'); hold on; plot(tiPred,zpred(1,:))


%% comparision between updated model and direct solution
ufun = @(t) 0.1*exp(-t)-1/(t+2);
FForced = @(t,x) F(t,x) + [zeros(n,1); (M\fExt)*ufun(t)];
opts = odeset('RelTol', 1e-4); 
[t,x] = ode15s(@(t,x) FForced(t,x),...
    linspace(0,2,1000), IC,opts);
xmid = (x(:,1:n)*phimid).';

% Check dynamics
% x0 = IMInfo.chart.map(IC);
x0 = [0;0];
[tiPred,xiPred] = ode45(@(t,x) R(t,x,ufun(t)),t,x0); 
tiPred = transpose(tiPred); xiPred = transpose(xiPred); 
% figure; plot(tiPred,xiPred)
paramFun = IMInfo.parametrization.map;
zpred = paramFun(xiPred);

figure; plot(t,xmid,'r-'); hold on; plot(tiPred,zpred(1,:))

%%
RDLin = RDInfo.reducedDynamics.coefficients(:,1:2);
Vmat  = IMInfo.parametrization.tangentSpaceAtOrigin;
Tmat  = Br*Vmat(1,:);

Tmp = RDInfo.eigenvectorsLinPart;
That = Tmp\(Tmat*Tmp);
LAMD = RDInfo.eigenvaluesLinPartFlow;

(That(1,1)+That(2,2))/(LAMD(1)+LAMD(2))

%% control stabilization via spectrum analysis of linearized system
% scheme 1
detR = det(RDLin);
b1   = Br(1);
b2   = Br(2);
r11  = RDLin(1,1); r12 = RDLin(1,2);
r21  = RDLin(2,1); r22 = RDLin(2,2);
pfun = @(k1,k2) r11+r22+b1*k1+b2*k2;
qfun = @(k1,k2) detR+(b1*r22-b2*r12)*k1+(b2*r11-b1*r21)*k2;

figure; hold on;
fimplicit(pfun,'--g','LineWidth',1.5);
fimplicit(qfun,[-100 100 -100 100],'-b','LineWidth',1.5);
set(gca,'FontSize',14);
xlabel('$k_1$','Interpreter','latex');
ylabel('$k_2$','Interpreter','latex');
% we infer that the third quadrant is feasible region
% so we can take k1=-100 and k2=-100
k1 = -100;
k2 = -100;
ks = [k1; k2];
eig(RDLin+Br*ks')

%% apply control law to stablize the unstable fixed pts
% comparision between updated model and direct solution
ufun = @(x) ks'*(Vmat.'*(repmat(phimid.'*x,[5,1])));
Fcontrol = @(t,x) F(t,x) + [zeros(n,1); (M\fExt)*ufun(x(1:nmodes))];
opts = odeset('RelTol', 1e-4); 
ICp = IC+[rand(nmodes,1)*0.05; zeros(n,1)]; % perturb all four modes
% ICp = IC; ICp(1,1) = rand(1)*0.05;
tf  = 30;
[tc,xc] = ode15s(@(t,x) Fcontrol(t,x),...
    linspace(0,tf,tf*1000), ICp,opts);
% without control
[tw,xw] = ode15s(@(t,x) F(t,x), linspace(0,tf,tf*1000), ICp,opts);
xcmid = (xc(:,1:n)*phimid).';
xwmid = (xw(:,1:n)*phimid).';
figure; plot(tc,xcmid,'r-'); hold on; plot(tw,xwmid,'b-');
legend('controlled','no control')

% Check dynamics
% x0 = IMInfo.chart.map(ICp);
x0 = Vmat.'*(repmat(phimid.'*ICp(1:nmodes),[5,1]));
ufeed = @(x) ks.'*x;
[tiPred,xiPred] = ode45(@(t,x) R(t,x,ufeed(x)),tc,x0); 
tiPred = transpose(tiPred); xiPred = transpose(xiPred); 
% figure; plot(tiPred,xiPred)
paramFun = IMInfo.parametrization.map;
zpred = paramFun(xiPred);
hold on; plot(tiPred,zpred(1,:))


%% scheme 2
tmp1 = That(1,1)+That(2,2)
tmp2 = That(2,2)-That(1,1)
% we observe that both tmp1 and tmp2 are negative. Thus, we are not able to
% find such control gain to stablize the unstable origin.
cgain = -0.01;
eig(diag(LAMD)+cgain*That)