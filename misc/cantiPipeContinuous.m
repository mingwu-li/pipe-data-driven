%-----------------------------------%
% BEGIN: hyperSensitiveContinuous.m %
%-----------------------------------%
function phaseout = cantiPipeContinuous(input)

mat = input.auxdata.mat;
Q  = mat.Q;
R  = mat.R;
C  = mat.C;
Br = mat.Br;
Rauto = input.auxdata.Rauto;
Wmap  = input.auxdata.Wmap;
t = input.phase.time;
x = input.phase.state;
u = input.phase.control;

% to adapt for the input structure of Rauto and Wmap, we first tranpose x
% and u
x = x.';
u = u.';
xdot = Rauto(x)+Br*u;
y = Wmap(x); 
z = C*y; 
run_cost = sum(z.*(Q*z),1)+sum(u.*(R*u),1);
phaseout.dynamics = xdot.';
phaseout.integrand = run_cost.';

%---------------------------------%
% END: hyperSensitiveContinuous.m %
%---------------------------------%
