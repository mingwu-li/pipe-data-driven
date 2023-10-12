function FRC = cont_SN_FRC(IMInfo,RDInfo,epsfrange,Omrange,run,lab,oid,amplitudeFunction,ispolar,varargin)


%% continuation of equilibria of reduced dynamics
% continuation in f such that it reaches the desired period
% call ep-toolbox
prob = coco_prob();
if ~isempty(varargin)
    prob = coco_set(prob, 'cont', 'h_max', varargin{1});
    prob = coco_set(prob, 'cont', 'NSV', 5);
end
prob = ode_ep2SN(prob, '', coco_get_id(run, 'ep'), lab);
% define monitor functions to state variables
[prob, args1, args2] = monitor_states(prob, ispolar);
runid = coco_get_id(oid, 'ep');
fprintf('\n Run=''%s'': Continue equilibria along primary branch.\n', ...
  runid);
coco(prob, runid, [], 1, {'om','eps',args1,args2}, {Omrange,epsfrange});


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

%% map back to physical coordinates
nlab = numel(om);
T    = RDInfo.transformation.map;
SSMFunction = IMInfo.parametrization.map;
phi = linspace(0,2*pi,128);
Aout = zeros(nlab,1);
%reduced coordinate
for j=1:nlab
    state = redSamp.z(j,:);
    state = [state conj(state)];
    zEval = transpose(state).*exp(1i*(phi));
    eta   = T(zEval); Zout = SSMFunction(eta);
    Zout  = Zout(3,:);
    Aout(j) = norm(Zout,'inf');
end
FRC = redSamp;
FRC.Aout = Aout;


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

