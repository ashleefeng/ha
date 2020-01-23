data1 = fret1;
data2 = fret2;

x1 = mean(data1);
x2 = mean(data2);
n1 = length(data1);
n2 = length(data2);

s1 = std(data1);
s2 = std(data2);
s = sqrt(s1^2/n1 + s2^2/n2);

temp1 = s1^2/n1;
temp2 = s2^2/n2;

df = (temp1 + temp2)^2/(temp1^2/(n1-1) + temp2^2/(n2-1));
fprintf('df: %f\n', df);

t = (x1 - x2)/s;
fprintf('t: %f\n', t);
