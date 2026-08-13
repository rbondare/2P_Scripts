function violin_kde(ax, xc, data, halfwidth, faceColor)
% VIOLIN_KDE  Toolbox-free violin: manual Gaussian KDE (Silverman bandwidth)
% mirrored about xc, with an IQR bar and a median dot. Degrades to a strip
% for tiny n. Shared by compare_plots and group_plots (one copy on purpose --
% see the heatmap_sorted duplication note in PIPELINE.md).
    data = data(isfinite(data));
    if numel(data) < 2
        if ~isempty(data); plot(ax, xc, data, 'o', 'Color', faceColor); end
        return;
    end
    s = std(data); if s <= 0; s = eps; end
    bw = 1.06 * s * numel(data)^(-1/5);            % Silverman's rule
    yq = linspace(min(data), max(data), 100);
    f = zeros(size(yq));
    for j = 1:numel(yq)
        f(j) = sum(exp(-0.5*((yq(j)-data)/bw).^2));
    end
    f = f / max(f) * halfwidth;                    % scale to halfwidth
    fill(ax, [xc + f, fliplr(xc - f)], [yq, fliplr(yq)], faceColor, ...
        'FaceAlpha', 0.5, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.5);
    q = prctile(data, [25 50 75]);
    plot(ax, [xc xc], [q(1) q(3)], 'k-', 'LineWidth', 2);
    plot(ax, xc, q(2), 'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
end
