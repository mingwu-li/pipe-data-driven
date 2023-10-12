function idx = extrema_index(x)

dx = x(2:end)-x(1:end-1);
dx = dx(2:end).*dx(1:end-1);
idx = find(dx<0);

end