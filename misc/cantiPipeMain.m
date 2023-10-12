function output = cantiPipeMain(Mat,Rauto,Wmap,eta0,tf)

% Parameters:
%-------------------------------------------------------------------%
%-------------------- Data Required by Problem ---------------------%
%-------------------------------------------------------------------%
auxdata       = struct();
auxdata.mat   = Mat;
auxdata.Rauto = Rauto;
auxdata.Wmap  = Wmap;

%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%
t0 = 0;
x0 = eta0; 
xMin = -50; 
xMax = +50; 
uMin = -100; 
uMax = +100; 

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0; 
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf; 
bounds.phase.finaltime.upper = tf;
bounds.phase.initialstate.lower = x0.'; 
bounds.phase.initialstate.upper = x0.';
bounds.phase.state.lower = [xMin,xMin]; 
bounds.phase.state.upper = [xMax,xMax];
bounds.phase.finalstate.lower = [xMin,xMin]; 
bounds.phase.finalstate.upper = [xMax,xMax]; 
bounds.phase.control.lower = uMin; 
bounds.phase.control.upper = uMax;
bounds.phase.integral.lower = 0; 
bounds.phase.integral.upper = 10000;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time    = [t0; tf]; 
guess.phase.state   = [x0.'; zeros(1,2)];
guess.phase.control = [0; 0];
guess.phase.integral = 0;

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name = 'canti-pipe-optimal-control';
setup.functions.continuous = @cantiPipeContinuous;
setup.functions.endpoint = @cantiPipeEndpoint;
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver = 'snopt';
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
% setup.scales.method = 'automatic-bounds';
setup.mesh.method = 'hp1';
setup.mesh.tolerance = 1e-6;
setup.mesh.maxiteration = 45;
setup.mesh.colpointsmax = 4;
setup.mesh.colpointsmin = 10;
setup.mesh.phase.colpoints = 4*ones(1,10);
setup.mesh.phase.fraction =  0.1*ones(1,10);
setup.method = 'RPMintegration';

%-------------------------------------------------------------------------%
%------------------------- Solve Problem Using GPOP2 ---------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);
solution = output.result.solution;

%--------------------------------------------------------------------------%
%------------------------------- Plot Solution ----------------------------%
%--------------------------------------------------------------------------%
figure;
pp = plot(solution.phase.time,solution.phase.state,'-o');
xl = xlabel('time');
yl = ylabel('state');
set(pp,'LineWidth',1.25,'MarkerSize',8);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
grid on

figure;
pp = plot(solution.phase.time,solution.phase.control,'-o');
xl = xlabel('time');
yl = ylabel('control');
set(pp,'LineWidth',1.25,'MarkerSize',8);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);

yout = Wmap(transpose(solution.phase.state));
zout = Mat.C*yout;
figure;
plot(solution.phase.time,zout)

end
