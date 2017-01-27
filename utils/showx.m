function showx( X, thresh, ncmap )
%SHOWX  Shows activation map values larger than threshold.

    keep = abs(X) > thresh*max(abs(X(:)));
    [keepi, keepj] = find(keep);
    colors = jet(ncmap);
    Xmax = max(max(X(:)), 0); Xmin = min(min(X(:)), 0);
    delta = ( Xmax - Xmin ) / (ncmap-1);
    
    nkeep = numel(keepi);
    cidx = fix( ( X(keep) - Xmin ) / delta ) + 1;

    hold on;
    for idx = 1:nkeep
        putcircle(keepj(idx), 255-keepi(idx), colors(cidx(idx),:));
    end
    hold off;
    colormap(colors); colorbar; caxis([Xmin Xmax]);
    xlim([1 size(X,2)]); ylim([1 size(X,1)]);
end

function putcircle( i, j, color )
    h = plot(i, j, 'o', 'Color', color);
    set(h,  'MarkerEdgeColor','k', ...
            'LineWidth', 1, ...
            'MarkerFaceColor',color);
end
