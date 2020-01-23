n1 = 66;
s1 = 0.1;
x1 = 0.5;
n2 = 32;
s2 = s1;
x2 = x1;
fret1sim = zeros(1, n1);
fret2sim = zeros(1, n2);
for i = 1:n1
    fret1sim(i) = normrnd(x1, s1);
end

for i = 1:n2
    fret2sim(i) = normrnd(x2, s2);
end

figure
histogram(fret1sim);
hold on
histogram(fret2sim);