function [M, C, K, fnl, fext, fphi, phinode] = buildModel(n,fc,u,type)

% set_parameters;
alpha = 0.001;
sigma = 0.000;
beta  = 0.5;
Acal  = 1000;
% u = 1;

% formulate mass, damp, gyro and stiff matrix as well as cubic nonlinear
% vector
mass  = eye(n);
damp  = zeros(n);
gyro  = zeros(n);
stiff = zeros(n);
fext  = zeros(n,1);
fphi  = zeros(n,1);
phinode = zeros(n,1);
cubic_coeff = zeros(n,n,n,n);
for i=1:n
    damp(i,i)  = sigma+alpha*(i*pi)^4;
    stiff(i,i) = (i*pi)^4-u^2*(i*pi)^2;
    switch type
        case 'first'
            fext(i) = fc*sqrt(2)*sin(i*pi*0.5);
            phinode(i) = sqrt(2)*sin(i*pi*0.5); % shape function at mid point
        case 'second'
            fext(i) = fc*sqrt(2)*(sin(i*pi*0.25)-sin(i*pi*0.75));
            phinode(i) = sqrt(2)*sin(i*pi*0.25); % shape function at 1/4 point
        otherwise
            error('type should be first or second');
    end
    fphi(i) = sqrt(2)*(1-(-1)^i)/(i*pi);
    for j=1:n
        if j~=i
            gyro(i,j) = 4*sqrt(beta)*u*i*j/(i^2-j^2)*(1-(-1)^(i+j));
        end
        cubic_coeff(i,j,j,i) = i^2*j^2;
    end
    
end

disp('the first four eigenvalues for undamped system');
lamd = eigs([zeros(n) eye(n);-mass\stiff -mass\(gyro)],4,'smallestabs')

% viscoelastic damping with nonlinear part
tmp   = sptensor(cubic_coeff);
subs1 = tmp.subs;
subs2 = subs1;
subs2(:,4) = subs1(:,4)+n;
subs = [subs1; subs2];
vals = [0.5*Acal*tmp.vals; alpha*Acal*tmp.vals];
f3 = pi^4*sptensor(subs, vals, [n,2*n,2*n,2*n]); 

M = mass;
C = damp+gyro;
K = stiff;
fnl = {[],f3};


end