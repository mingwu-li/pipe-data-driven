function fred = calibrateQFRC(RDInfo,Ts,Omega)
% CALIBRATEQFRC This function calibrate the forcing in reduced dynamics via
% the second frequency component of a quasi-periodic signal

coeffs = RDInfo.conjugateDynamics.coefficients;
alphas = real(coeffs);
omegas = imag(coeffs);
% locate root of alpha(rho)
polycoeffs = zeros(1,2*numel(alphas)-1);
polycoeffs(1:2:end) = alphas(end:-1:1);
polycoeffs(2:2:end-1) = 0;
rhos = roots(polycoeffs);
rhos = sort(rhos);
idx  = find(imag(rhos)==0);
assert(~isempty(idx),'no limit cycle is found');
rhos = abs(rhos(idx(1)));
omgs = 0;
for k=1:numel(omegas)
    omgs = omgs+omegas(k)*rhos^(2*k-2);
end
fred = zeros(numel(Ts),1);
% loop over Ts (Omega)
for k=1:numel(Ts)
    omp  = omgs-Omega(k);
    t0 = linspace(0,2*pi/omp,128)';
    x0 = [rhos*cos(omp*t0) rhos*sin(omp*t0)];
    p0 = [0,Omega(k)];
    % continuation in f such that it reaches the desired period
    odefunc = @(x,p) conjDynCart(x,p,coeffs,1);
    prob = coco_prob();
    prob = coco_set(prob, 'cont', 'NAdapt', 2);
    prob = ode_isol2po(prob, '', odefunc, t0, x0, {'f','Om'}, p0);
    runk = ['cal_f',num2str(k)];
    bd = coco(prob, runk, [], 1, {'po.period','f','Om'}, {[0,Ts(k)],[0,1000]});
    % record corresponding f
    eplab = coco_bd_labs(bd,'EP');
    eplab = max(eplab);
    sol = po_read_solution('', runk, eplab);
    fred(k) = sol.p(1);
end

end

