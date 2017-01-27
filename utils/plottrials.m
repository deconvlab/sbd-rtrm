load('trials.mat'); clf;


ha = tight_subplot(3, trials+1, [.02 .02], [.01 .05]+.05, [.03 .02]+.05);
plottrial(ha, trials, 0, A0, X0, Y, []);
for i = 1:trials
    plottrial(ha, trials, i, Aouts{i}, Xouts{i}, cconvfft2(Aouts{i}, Xouts{i}), lambda(i));
end