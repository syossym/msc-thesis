function [min_mat, max_mat] = GetExtrema(x,y,num,h_Extrema)

% [detected_mins, detected_maxs] = peakdet(y, delta, x);
y = smooth(y);

[max_indices_x, max_indices_y] = GetMaximumIndices(x,y);
detected_maxs = [max_indices_x.', max_indices_y.'];

V = 1-y;
[min_indices_x, min_indices_y, indices]  = GetMaximumIndices(x,V);
detected_mins = [min_indices_x.', min_indices_y.'];

[s_maxs, i_maxs] = sort(detected_maxs(:,2), 1, 'descend');
[s_mins, i_mins] = sort(detected_mins(:,2), 1, 'descend');

sorted_mins = detected_mins(i_mins, :);
sorted_maxs = detected_maxs(i_maxs, :);
min_mat = sorted_mins(1:min([num,length(sorted_mins(:,1))]), :);
max_mat = sorted_maxs(1:min([num,length(sorted_maxs(:,1))]), :);
min_mat(:,2) = 1 - min_mat(:,2);

[s_maxs, i_maxs] = sort(max_mat(:,1), 1, 'ascend');
[s_mins, i_mins] = sort(min_mat(:,1), 1, 'ascend');
max_mat = max_mat(i_maxs, :);
min_mat = min_mat(i_mins, :);

% Plotting
figure(h_Extrema); box on;
plot(x,y,'b'); hold on;
for (mn=1:length(max_mat(:,1)))
    plot(max_mat(mn,1), max_mat(mn,2), 'rx')
    text(max_mat(mn,1), max_mat(mn,2), num2str(mn));
end
for (mn=1:length(min_mat(:,1)))
    plot(min_mat(mn,1), min_mat(mn,2), 'r*');
    text(min_mat(mn,1), min_mat(mn,2), num2str(mn));
end
hold off;
drawnow;

function [M_x,M_y,indices] = GetMaximumIndices(x,V)

A = V;

interval = length(A)/40;
d_conv = [1 -1];                         % derivative kernel
dif = conv(A, d_conv);                   % first derivative
dif = dif(1:(length(dif)-2));
dif_dif = conv(dif, d_conv);             % second derivative
dif_dif = dif_dif(1:(length(dif_dif)-2));

a = [];
for ii = 2:length(A)-1
    if ((A(ii)-A(ii-1))>0 & (A(ii+1)-A(ii))<0)
        a = [a ii];
    end
end

indices = [];
M = zeros(1, length(V));
for ii=1:length(M)
    for jj=1:length(a)
        if ii==a(jj)
            M_y(ii)=V(a(jj));
            M_x(ii)=x(a(jj));
            indices = [indices, a(jj)];
        end
    end
end

M_y = M_y(M_y~=0);
M_x = M_x(M_x~=0);

