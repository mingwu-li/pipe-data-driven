function  y = conjDynPolar(x,p,coeffs,kapa)
% CONJDYNCART This function returns the vector field of conjugate reduced
% dynamics in polar coordinates.

alphas = real(coeffs);
omegas = imag(coeffs);
rho = x(1,:);
psi = x(2,:);
epf = p(1,:);
Om  = p(2,:);
y   = zeros(2,numel(rho));
as  = zeros(1,numel(rho));
os  = zeros(1,numel(rho));
for k=1:numel(alphas)
    as = as+alphas(k)*rho.^(2*k-2);
    os = os+omegas(k)*rho.^(2*k-2);
end
y(1,:) = rho.*as+kapa*Om.^2.*epf.*sin(psi);
y(2,:) = os-Om+kapa*Om.^2.*epf.*cos(psi)./rho;

% y(1,:) = rho.*as+f.*sin(psi);
% y(2,:) = os-Om+f.*cos(psi)./rho;

end