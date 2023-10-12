function FRC = cont_po_FRC(IMInfo, RDInfo, kapa, epsf, Omrange,oid, amplitudeFunction, ispolar, varargin)


%% continuation of equilibria of reduced dynamics
coeffs = RDInfo.conjugateDynamics.coefficients;
% continuation in f such that it reaches the desired period
if ispolar
    odefun = @(z,p) conjDynPolar(z,p,coeffs,kapa);
else
    odefun = @(z,p) conjDynCart(z,p,coeffs,kapa);
end
prob = coco_prob();
prob = coco_set(prob, 'corr', 'ItMX', 20);
% prob = coco_set(prob, 'cont', 'PtMX', [80,0]);
if ~isempty(epsf)
    p0   = [epsf; Omrange(1)];
    [p0,z0] = initial_fixed_point(p0,ispolar,odefun);
end
if ~isempty(varargin)
    prob = coco_set(prob, 'cont', 'h_max', varargin{1});
    if numel(varargin)>1 
        p0 = varargin{2}; 
        [p0,z0] = initial_fixed_point(p0,ispolar,odefun,varargin{3});
    end
end

% call ep-toolbox
prob = ode_isol2ep(prob, '', odefun, z0, {'eps','om'}, p0);
% define monitor functions to state variables
[prob, args1, args2] = monitor_states(prob, ispolar);
runid = coco_get_id(oid, 'ep');
fprintf('\n Run=''%s'': Continue equilibria along primary branch.\n', ...
  runid);
coco(prob, runid, [], 1, {'om',args1,args2,'fred'}, Omrange);


%% extract results
bd = coco_bd_read(runid);
om  = coco_bd_col(bd, 'om')';
if ispolar
    rho = coco_bd_col(bd, args1)'; 
    th  = coco_bd_col(bd, args2)';
    state = rho.*exp(1j*th);
else
    zRe = coco_bd_col(bd, args1)'; 
    zIm = coco_bd_col(bd, args2)';
    rho = sqrt(zRe.^2+zIm.^2);
    th  = atan2(zIm,zRe);
    state = zRe+1j*zIm;
end
epsf = coco_bd_col(bd, 'eps')';
stab = coco_bd_col(bd, 'eigs')'; 
stab = real(stab);
stab = all(stab<0,2);
SNidx = coco_bd_idxs(bd, 'SN');
HBidx = coco_bd_idxs(bd, 'HB');
% covert th such that it is in [0,2pi] and rho such it is postive
th             = mod(th,2*pi);
negRhoIdx      = find(rho<0);
rho(negRhoIdx) = -rho(negRhoIdx);
th(negRhoIdx)  = th(negRhoIdx)-pi;

redSamp     = struct();
redSamp.om  = om;
redSamp.z   = state;
redSamp.rho = rho;
redSamp.th  = th;
redSamp.ep  = epsf;
redSamp.st  = stab;
redSamp.SNidx = SNidx;
redSamp.HBidx = HBidx;

%% map back to physical coordinates
nlab = numel(om);
T    = RDInfo.transformation.map;
SSMFunction = IMInfo.parametrization.map;
phi = linspace(0,2*pi,128);
Aout = zeros(nlab,1);
%reduced coordinate
for j=1:nlab
    state = redSamp.z(j,:);
    zEval = state*exp(1i*(phi));
    eta   = T([zEval; conj(zEval)]); 
    Zout  = SSMFunction(eta);
    Zout  = Zout(3,:);
    Aout(j) = norm(Zout,'inf');
end
FRC = redSamp;
FRC.Aout = Aout;


end


function [p0,z0] = initial_fixed_point(p0,ispolar,odefun,varargin)
% INITIAL_FIXED_POINT This function construct initial solution to the fixed
% point of leading-order reduced dynamics. Two methods: forward simulation
% and optimization, are avaliable to obtain such an initial fixed point.

if isempty(varargin)
    if ispolar
        z0 = 0.1*ones(2,1);
    else
        z0 = zeros(2,1);
    end
else
    z0 = varargin{1};
end
% construct initial guess equilibrium points
fsolveOptions = optimoptions('fsolve','MaxFunctionEvaluations',1000000,...
    'MaxIterations',1000000,'FunctionTolerance', 1e-10,...
    'StepTolerance', 1e-8, 'OptimalityTolerance', 1e-10);
z0 = fsolve(@(z) odefun(z,p0),z0,fsolveOptions);

if ispolar % regularize initial solution if it is in polar form
    z0(2:2:end) = mod(z0(2:2:end),2*pi); % phase angles in [0,2pi]
    if z0(1)<0
        z0(1) = -z0(1);      % positive amplitudes
        z0(2) = z0(2)+pi;
    end
end
end

function [prob, args1, args2] = monitor_states(prob, ispolar)
% MONITOR_STATES This function add state of reduced dynamics as
% continuation parameters and define string array for the state. The prob
% here is a continuation problem. (args1,args2)=(rhoargs,thargs) or 
% (Reargs,Imargs) depending on the value of ispolar

if ispolar
    args1 = 'rho';
    args2 = 'psi';
    prob = coco_add_pars(prob, 'radius', 1, args1);
    prob = coco_add_pars(prob, 'angle', 2, args2);
else
    args1 = 'Rez';
    args2 = 'Imz';
    prob = coco_add_pars(prob, 'realPart', 1, args1);
    prob = coco_add_pars(prob, 'imagPart', 2, args2);
end 

end

