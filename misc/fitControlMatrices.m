function [Br,regErrors] = fitControlMatrices(yuxData,RDInfo,regBr)
% Assume graph type approach where reduced coordinates are a linear
% projection. xuTraj is a list of full trajectories (of potentially delay 
% embedded data) where first column is time, second are the embedded
% observables, third the control and last the reduced coordinates 
% This code is adapted from the same function at 
% https://github.com/StanfordASL/SSMR-for-control 

redynFun = RDInfo.reducedDynamics.map;
t = [];       % time values
Xr = [];      % reduced coordinates at time k
U = [];       % controls at time k
dXrdt = [];   % time derivatives of reduced coordinates at time k
apprxOrd = 3; % approximation order for the derivative
% Data in matrices
kidx = 0;
for ii = 1:size(yuxData,1)
    t_in = yuxData{ii,1}; Y_in = yuxData{ii,2}; 
    U_in = yuxData{ii,3}; Xr_in = yuxData{ii,4}; 
    [dXridt,Xri,ti] = finiteTimeDifference(Xr_in,t_in,apprxOrd);
    Ui = U_in(:, 1+apprxOrd:end-apprxOrd);
    t = [t ti]; Xr = [Xr Xri]; dXrdt = [dXrdt dXridt]; 
    U = [U Ui]; kidx = [kidx numel(t)];
end

% Fit reduced dynamics control matrix
deltaDerivatives = dXrdt - redynFun(Xr);
% Learn whole B matrix
[Br,~,~] = ridgeRegression(U, deltaDerivatives, ones(size(t)), [], regBr);
regErrorBr = mean(sqrt(sum((deltaDerivatives - Br*U).^2)))/...
             max(sqrt(sum((deltaDerivatives).^2)));
         
% Output errors of the regression in %
regErrors = regErrorBr*100;

% Plot reconstructed rhs and dXrdt
gap = max(1, floor(size(yuxData,1)/4));
for k=1:gap:size(yuxData,1)-1
    figure;
    idx = kidx(k)+1:kidx(k+1);
    plot(t(idx),dXrdt(1,idx)); hold on
    plot(t(idx),dXrdt(2,idx)); 
    rhs = redynFun(Xr)+Br*U;
    plot(t(idx),rhs(1,idx),'r--'); hold on
    plot(t(idx),rhs(2,idx),'m--'); 
    % Integration of xi
    R = @(t,x,u) redynFun(x) + Br*u;
    ufun = @(ts) interp1(t(idx),U(idx),ts);
    [tiPred,xiPred] = ode45(@(t,x) R(t,x,ufun(t)),t(idx),Xr(:,idx(1))); 
    tiPred = transpose(tiPred); xiPred = transpose(xiPred); 
    figure; plot(tiPred,xiPred); hold on; 
    plot(t(idx),Xr(:,idx))
end




end