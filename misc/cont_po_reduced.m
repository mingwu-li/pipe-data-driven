function runom = cont_po_reduced(RDInfo,kapa,fFull,omegaSpan,runid,varargin)
% CONT_PO_REDUCED This function perform continuation of limit cycles of the
% reduced dynamics. We first perform continuation in forcing amplitude to
% reach the desired forcing level. We then perform continuation in forcing
% frequency to yield a family of periodic orbits. Since we have two
% families of periodic orbits along the left and right side. The above
% procedure will be peformed two times.

% setup
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
% loop over left and right
omend = [omegaSpan(1) omegaSpan(end)];
runf  = {'contFleft','contFright'};
runom = {['contOmleft_',runid],['contOmright_',runid]};
for k=1:2
    omp  = omgs-omend(k);
    t0 = linspace(0,2*pi/abs(omp),128)';
    x0 = [rhos*cos(omp*t0) rhos*sin(omp*t0)];
    p0 = [0,omend(k)];
    % continuation in forcing amplitude such that it reaches the desired
    % forcing level
    odefunc = @(x,p) conjDynCart(x,p,coeffs,kapa);
    prob = coco_prob();
    prob = coco_set(prob, 'cont', 'NAdapt', 2, 'NSV', 5);
    prob = coco_set(prob, 'corr', 'ItMX', 20);
    prob = ode_isol2po(prob, '', odefunc, t0, x0, {'f','Om'}, p0);
    bd1  = coco(prob, runf{k}, [], 1, {'f','po.period','Om'}, [0,fFull]);
    % check whether the desired forcing level is reached
    eplab = coco_bd_labs(bd1,'EP');
    eplab = max(eplab);
    sol   = po_read_solution('', runf{k}, eplab);
    assert(abs(sol.p(1)-fFull)<1e-3*fFull,'desired forcing level is not reached');
    % continuation in forcing frequency to yield a family of periodic
    % orbits
    prob = coco_prob();
    prob = coco_set(prob, 'cont', 'NAdapt', 2, 'NSV', 5, 'h_max',5);
    prob = ode_po2po(prob, '', runf{k}, eplab);
    if ~isempty(varargin)
        prob = coco_add_event(prob, 'UZ', 'Om', varargin{1}{k});
    end
    bd2  = coco(prob, runom{k}, [], 1, {'Om','po.period','f'}, omegaSpan);
    figure; hold on
    coco_plot_bd(runom{k},'Om','po.period');
end

end

