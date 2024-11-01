function PlotWaterfall(x, y_mat, handler, color, x_grid_range, y_var, coeff, scale)

if (nargin < 8)
    scale = 1;
end

[h,w] = size(y_mat);
loc = [0:1:h-1];

if (nargin < 7)
    coeff = 1;
end

try
    axes(handler);
catch
    figure(handler);
end

scale = 2;

for (ii=1:h)
    hold on; box on;
    if (scale==1)
        y_norm = coeff*y_mat(ii,:)/max(y_mat(ii,:));
    elseif(scale == 2)
        y_norm = coeff*y_mat(ii,:)./max(max(y_mat));
    else
        y_norm = coeff*y_mat(ii,:);
    end
    y = loc(ii) + 3.*y_norm;
    if (mean(y_norm) > 0.5)
        y = y - mean(y_norm);
    end
    plot(x, y, 'Color', color, 'LineWidth', 0.2);
    scaling = max(max(y_mat))./max(coeff*y_mat(ii,:));
    scaling_text = ['{\color{' color '}\times' num2str(round(scaling)) '}'];
    [m,n] = max(y); 
    %text(x(n),y(n),scaling_text);
    if (nargin >= 6)
        if (~isnan(y_var.value))
            if (strcmp(y_var.num_type, 'exp'))
                text(x_grid_range(1), y(1) + 0.2, [y_var.name '=' strrep(num2str(y_var.value(ii),'%1.1e'), 'e+0', 'x10^{') '}' y_var.units]);
            else
                text(x_grid_range(1), y(1) + 0.2, [y_var.name '=' num2str(y_var.value(ii)) y_var.units]);
            end
        end
    end
    %text(x_grid_range(end)-0.01, y(end) + 0.2, [', \times' num2str(scaling, '%1.2f')]);
end
axis tight;
set(gca, 'XLim', x_grid_range);
set(gca, 'YTickLabel', '');