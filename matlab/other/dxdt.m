function [ f ] = dxdt( t,a )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

A1=a(1);
A2=a(2);
A3=a(3);
A4=a(4);
%A4=a(4);
%ZZ=x(5);

p=[.1,.03,.3,.01];

dAAdt = -p(1)*A1-p(2)*A1;
dAZdt = p(1)*A1-p(3)*A2;
dZAdt = p(2)*A1-p(4)*A3;
dZZdt = p(3)*A2 + p(4)*A3;

f = [dAAdt; dAZdt; dZAdt; dZZdt];

end

