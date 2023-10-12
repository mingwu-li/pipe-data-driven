function  FRC = quasi_po_FRC(IMInfo,RDInfo,runid,labs)
% QUASI_PO_FRC This function maps limit cycles in reduced system to tori in
% the full system via a coordinate transformation. For each limit cycle, we
% first construct the corresponding torus in reduced coordinates and then
% map it to physical coordinates.


% setup
nt   = 128;
nlab = numel(labs);
sols = cell(nlab,1);
stab = false(nlab,1);
omeg = zeros(nlab,1);
epsf = zeros(nlab,1);
for k = 1:nlab
    sol = po_read_solution('',runid, labs(k));
    stab(k) = all(abs(sol.po_test.la)<1);
    omeg(k) = sol.p(2);
    epsf(k) = sol.p(1);
    sols{k} = sol;
end

FRC = struct();
FRC.om  = omeg;
FRC.st  = stab;
FRC.ep  = epsf;
FRC.lab = labs;

%% construct torus in reduced system using interpolation of periodic cycle
disp('Constructing torus in reduced dynamical system');
qTr  = cell(nlab,1);
tTr  = cell(nlab,1);
nSeg = zeros(nlab,1);
dim  = size(sol.xbp,2);
for i=1:nlab
    soli = sols{i};
    tbp = soli.tbp;
    xbp = soli.xbp;
    om  = soli.p(2);
    Ti  = soli.T;
    assert(abs(om-omeg(i))<1e-3*omeg(i), 'Read wrong solution from data');
    fprintf('Interpolation at (omega,epsilon) = (%d,%d)\n', om, soli.p(2));
    tsamp   = linspace(0,2*pi/om,nt);
    numSegs = numel(tbp);
    ys = zeros(nt,dim,numSegs);
    for k=1:numSegs
        % shift tsamp such that the initial time is the same as
        % tbp(k)
        tsampk = tsamp+tbp(k);
        tsampk = mod(tsampk,Ti); % mod with the period of the periodic cycle        
        xk  = interp1(tbp, xbp, tsampk, 'pchip'); % Update basepoint values
        zk  = xk(:,1:2:end-1)+1j*xk(:,2:2:end);
        zk  = zk.*exp(1j*(om.*tsamp'));
        ys(:,1:2:end-1,k) = real(zk);
        ys(:,2:2:end,k)   = imag(zk);
    end
    qTr{i} = ys;
    tTr{i} = tsamp;
    nSeg(i) = numSegs;               
end      
FRC.qTr = qTr;
FRC.tTr = tTr;
FRC.nSeg = nSeg;

%% map it back to physical coordinates
zTr = cell(nlab,1);
dim = size(IMInfo.parametrization.tangentSpaceAtOrigin,1);
T   = RDInfo.transformation.map;
SSMFunction = IMInfo.parametrization.map;
Aout = zeros(nlab,1);
for j=1:nlab
    qTrj = qTr{j};
    tTrj = tTr{j};
    nt   = numel(tTrj);
    numSegs = nSeg(j);
    Zout_frc = zeros(nt,dim,numSegs);
    for k=1:numSegs
        xbp = qTrj(:,:,k);
        xbp = xbp';
        x_comp = xbp(1:2:end-1,:)+1i*xbp(2:2:end,:);
        zEval  = [x_comp; conj(x_comp)];
        eta    = T(zEval);
        Zout   = SSMFunction(eta);
        Zout_frc(:,:,k) = Zout';
    end
    zTr{j} = Zout_frc;
    tmp = Zout_frc(:,2,:);
    Aout(j) = norm(tmp(:),'inf');
end

FRC.zTr = zTr; 
FRC.Aout = Aout;
end