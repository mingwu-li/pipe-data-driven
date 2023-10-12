function  y = conjDynCart(x,p,coeffs,kapa)
% CONJDYNCART This function returns the vector field of conjugate reduced
% dynamics in cartesian coordinates.

alphas = real(coeffs);
omegas = imag(coeffs);
x1 = x(1,:);
x2 = x(2,:);
ep = p(1,:);
Om = p(2,:);
rho = x1.^2+x2.^2;
y   = zeros(2,numel(x1));
as  = zeros(1,numel(x1));
os  = zeros(1,numel(x1));
for k=1:numel(alphas)
    as = as+alphas(k)*rho.^(k-1);
    os = os+omegas(k)*rho.^(k-1);
end
y(1,:) = x1.*as-x2.*(os-Om);
y(2,:) = x2.*as+x1.*(os-Om)+kapa*Om.^2.*ep;

end