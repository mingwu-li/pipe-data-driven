%---------------------------------%
% BEGIN: hyperSensitiveEndpoint.m %
%---------------------------------%
function output = cantiPipeEndpoint(input)

mat = input.auxdata.mat;
Qf  = mat.Qf;
C   = mat.C;
Wmap = input.auxdata.Wmap;
x = input.phase.finalstate;
x = x.';
y = Wmap(x); 
z = C*y; 
q_may = z'*Qf*z; 
q_lag = input.phase.integral;
output.objective = q_may+q_lag;
%---------------------------------%
% BEGIN: hyperSensitiveEndpoint.m %
%---------------------------------%
