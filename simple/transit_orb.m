%% Finding a 2D SSM for a simply supported pipe conveying fluid: mixed mode case
% This is an example of how to reconstruct a slow 2D SSM of a mechanical system 
% using synthetic measurements of a scalar quantity. In this example, we measure 
% the midpoint displacement of a simply supported geometrically nonlinear pipe 
% conveying fluid. [1] In this livescript, we choose flowspeed=4 such that the 
% dynamics is in the supercritical region.
% 
% [1] Paidoussis, M. P. (1998). _Fluid-structure interactions: slender structures 
% and axial flow_ (Vol. 1). Academic press. 
%
% Here we take flow speed to be 4 to investigate post-buckling dynamics. In
% particular, we are interested in the dynamics on a mixed-mode SSM, where
% we can capture the coexistence of three static configurations: the
% straight, unstable configuration, and two symmetric buckled, stable
% configurations.

clearvars
close all
%% Example setup
nmodes = 4;
fcload = 1;
flowspeed = 4;
[M, C, K, fnl, fExt, fphi, phimid] = buildModel(nmodes,fcload,flowspeed,'first');
n = size(M,1);    % mechanical dofs (axial def, transverse def, angle)
[F, lambda] = functionFromTensors(M, C, K, fnl);
lambda = sort(lambda);
%% Generation of Synthetic Data
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex');
[V,D,W] = DS.linear_spectral_analysis();
% four ICs, where two of them go to buckled configuration directly, and the
% rest two undergo large-amplitude oscillation initially and then approach
% to the stable buckled static configurations.
loadAmp = [1e-4 -1e-4 1e1 -1e1 6e1 -6e1]; % amplitude of static forces
nTraj = numel(loadAmp);
IC = zeros(2*nmodes,nTraj);
for i = 1 : nTraj % calculate initial conditions
    IC(:,i) = getStaticResponse(K, M, F, fExt*loadAmp(i));
end    

newMeasurement = true;
observable = @(x) phimid'*x(1:n,:);
slowTimeScale = 2*pi/abs(lambda(1));
if newMeasurement
    numberPointsPerPeriod = 200; 
    endTime = 30;
    nSamp = endTime*numberPointsPerPeriod+1;
    dt = endTime/(nSamp-1);
    xData = integrateTrajectories(F, endTime, IC , nSamp, observable, 'odesolver', @ode45);
    save('transit.mat', 'xData', 'dt', 'endTime', 'nSamp', 'lambda')
else
    load transit.mat
end

%%
customFigure();
for k=1:nTraj
    plot(xData{k,1},xData{k,2}); hold on
end  
xlabel('$t \, [$s$]$','Interpreter','latex'); 
ylabel('$u \, [$m$]$','Interpreter','latex'); 
legend({'Trajectory 1', 'Trajectory 2'})
title('Generated data')

%% Delay embedding
SSMDim = 2;
overEmbed = 0;
yData = coordinatesEmbedding(xData, SSMDim, 'OverEmbedding', overEmbed);
%% Truncate signals that are far away from the SSM
indTrain = [5 6];
indTest  = [1 2 3 4];
sliceInt = [1*slowTimeScale, endTime]; yDataTrunc(1:2,:) = sliceTrajectories(yData(1:2,:), sliceInt);
sliceInt = [5*slowTimeScale, endTime]; yDataTrunc(3:6,:) = sliceTrajectories(yData(3:6,:), sliceInt);

%% Data-driven manifold fitting
SSMOrder = 1; % we fit a flat SSM here
[IMInfo, SSMChart, SSMFunction] = IMGeometry(yDataTrunc(indTrain,:), SSMDim, SSMOrder);

%%
etaData = projectTrajectories(IMInfo, yData);
etaDataTrunc = projectTrajectories(IMInfo, yDataTrunc);
%% 
% We plot the test and training set trajectories projected onto the plane. After 
% all initial transients die out, the trajectories seem to be confined to a slowly 
% decaying oscillating eigenmode. This is the data we will use to train our reduced 
% order model for the dynamics. In these coordinates, the initial transient that 
% we removed is visible as irregular high-frequency oscillations in the first 
% revolution.
plotReducedCoordinates(etaData);
legend({'Test set trajectory', 'Training set trajectory'})

plotReducedCoordinates(etaDataTrunc);
legend({'Test set trajectory-trunc', 'Training set trajectory-trunc'})

%% 
% Furthermore, we draw the end midpoint component of the manifold shape along with 
% the truncated trajectory from the training set. 
plotSSMWithTrajectories(IMInfo, 1, yDataTrunc(indTrain,:))
view(-100,20); zlabel('$u \, [$m$]$','Interpreter','latex')


%% Reduced order model
ROMOrder = 13; % we fit a reduced dynamics in polynomial form up to O(13)
RDInfo = IMDynamicsFlow(etaDataTrunc(indTrain,:), ...
    'R_PolyOrd', ROMOrder, 'style', 'modal');
%% 
[yRec, etaRec, zRec] = advect(IMInfo, RDInfo, yDataTrunc);
%% Evaluation of reduced dynamics
normedTrajDist = computeTrajectoryErrors(yRec, yDataTrunc);
NMTE = mean(normedTrajDist(indTest))

%% 
plotReducedCoordinates(etaDataTrunc(indTest,:), etaRec(indTest,:))
legend({'Test set (truncated)', 'Prediction'})
%% 
figure; 
plotTrajectories(yData(indTest(1),:), yRec(indTest(1),:), 'm', 'DisplayName', {'Test set', 'Prediction'})
ylabel('$u \, [$m$]$','Interpreter','latex')

%% 
plotTrajectories(yData(indTest([1,4]),:), yRec(indTest([1,4]),:), 'm', 'DisplayName', {'Test set', 'Prediction'})
ylabel('$u \, [$m$]$','Interpreter','latex')

