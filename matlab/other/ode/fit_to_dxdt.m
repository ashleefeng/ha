t = [0 1 5 10 20 45 60];
laa = [1 0.64358 0.45481 0.33923 0.22278 0.14511 0.15896];
laz = [0 0.18863 0.13806 0.10849 0.05921 0.00817 0.00668];
lza = [0 0.08501 0.12894 0.15909 0.16447 0.10388 0.08184];
lzz = [0 0.08278 0.27819 0.39319 0.55354 0.74284 0.75253];

x = zeros(1, 8);
x(1:4) = [.1,.03,.3,.01];
x(5:8) = [1,0,0,0];
lb = [0, 0, 0, 0, 0, 0, 0, 0];
ub = [1, 1, 1, 1, 1, 1, 1, 1];
ydata = [laa;laz;lza;lzz];
[pbest,presnorm,presidual,exitflag,output] ...
    = lsqcurvefit(@paramfun,x,t,ydata, lb, ub);

fprintf('New parameters: %f, %f, %f, %f\n',pbest(1:4));
fprintf('Original parameters: %f, %f, %f, %f\n',x(1:4));

figure
hold on
odesl = presidual + ydata;
plot(t, odesl, 'r');
plot(t, ydata, 'b');
%legend('ODE Solution','Experiment')
y=paramfun([pbest(1), pbest(2), pbest(3), pbest(4), 1, 0, 0, 0], [0 60]);
%plot(t,y, 'g');
hold off