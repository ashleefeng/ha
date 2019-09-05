files = dir('*.dat');
for file = files'
    trace = readtable(file.name);
    fret = trace.Var3 ./ (trace.Var2 + trace.Var3);
    h = histogram(fret);
    saveas(h, strcat(file.name, '.png'));
end
