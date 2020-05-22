function pos = paramfun(x, tspan)

p = x(1:4);
A1 = x(5);
A2 = x(6);
A3 = x(7);
A4 = x(8);

dAAdt = -p(1)*A1-p(2)*A1;
dAZdt = p(1)*A1-p(3)*A2;
dZAdt = p(2)*A1-p(4)*A3;
dZZdt = p(3)*A2 + p(4)*A3;

f = @(t, a) [dAAdt; dAZdt; dZAdt; dZZdt];

[~, pos] = ode15s(f, tspan, x(5:8));

