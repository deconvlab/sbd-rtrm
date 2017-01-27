function plottrial( ha, trials, trial, A, X, Y, lambda)
    plotatidx(ha, trials, trial, 1, A);
    if trial == 0
        ylabel('Kernel $$\mathcal{A}_0\ (\hat{\mathcal{A}})$$', ...
            'Interpreter', 'latex');
    end
    
    plotatidx(ha, trials, trial, 2, X);
    if trial == 0
        ylabel('Activation map $$\mathcal{X}_0\ (\hat{\mathcal{X}})$$', ...
            'Interpreter', 'latex');
    end
    
    plotatidx(ha, trials, trial, 0, Y);
    if trial == 0
        ylabel('Observation $$\mathcal{Y}$$ (reconstruction $$\hat{\mathcal{Y}})$$',...
            'Interpreter', 'latex');
        title('\bf Ground Truth', 'Interpreter', 'Latex');
    else
        title(sprintf('\\textbf{RTRM Estimate %d:} $\\lambda = %.03f$', trial, lambda), ...
            'Interpreter', 'latex');
    end
    
    drawnow;
end

function plotatidx( ha, trials, trial, idx, varin )
    h = ha(idx*(1+trials) + trial + 1);
    axes(h);
    imagesc(abs(varin));
    set(h, 'XTick', [], 'YTick', []);
end