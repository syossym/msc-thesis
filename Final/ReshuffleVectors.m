function F = ReshuffleVectors(E_min_MCQW, E_min_MC)

temp = E_min_MCQW;
temp(temp==0)=NaN;

[h,w]=size(temp);

flag = zeros(h,w);
new_temp = nan(h,w);

thr =6e-4;

for (jj=1:h)
    first_num_index = find(~isnan(temp(jj,:)),1,'first');
    %first_num_index = find(flag(jj,:)==0, 3, 'first');
    [i,j] = find(temp(:,first_num_index)~= new_temp(:,first_num_index));
    new_temp(jj,first_num_index)=temp(i(1),first_num_index);
    flag (jj,first_num_index)=1;
    kk=1;
    for (ii=first_num_index+1:w)
        
        %if (~isnan(new_temp(jj,ii-1)))
        dist_vec = abs(temp(:,ii)-new_temp(jj,ii-1));
        %else
        %    dist_vec = abs(temp(:,ii)-new_temp(jj,kk));
        %end
        
        if (isnan(min(dist_vec)))
            new_temp(jj,ii)=NaN;
            %flag (jj, ii)=1;
            continue;
        elseif (min(dist_vec)>=thr)
            new_temp(jj,ii)=NaN;
            kk=ii;
            %flag (jj, ii)=1;
            continue;
        end;
        
        [tt, index_min]=min(dist_vec);
        
        if (flag (index_min, ii)==0)
            new_temp(jj,ii)=temp(index_min, ii);
            flag (index_min, ii)=1;
        end;
    end;
end;

%new_temp(flag==0)=temp(flag==0)

%new_temp(new_temp==0)=NaN;

%new_temp = [new_temp;  resize(temp(flag==0), 1, length(new_temp(1,:)))];

flag(isnan(temp)) = nan;
new_temp(new_temp==0) = nan;

for (hh=1:h)
    if (length(flag(hh,flag(hh,:)==0))>20)
        vec = flag(hh, :);
        vec(vec==1)=nan;
        vec(vec==0)=temp(hh, flag(hh,:)==0);
        new_temp = [new_temp; vec];
    end
end

for (hh=1:length(new_temp(:,1)))
    nans(hh) = sum(isnan(new_temp(hh,:)));
end

new_temp = new_temp(nans<(3/4)*length(new_temp(1,:)), :);

for (hh=1:length(new_temp(:,1)))
    for (ww=2:w-1)
        if (~isnan(new_temp(hh,ww)) && isnan(new_temp(hh,ww-1)) && isnan(new_temp(hh,ww+1)))
            new_temp(hh,ww) = nan;
        end
    end
end

for (hh=1:length(new_temp(:,1)))
    first_non_nan = find(~isnan(new_temp(hh,:)), 1, 'first');
    last_non_nan = find(~isnan(new_temp(hh,:)), 1, 'last');
    new_temp(hh, 1:first_non_nan-1) =  new_temp(hh,first_non_nan);
    new_temp(hh, last_non_nan+1:end) =  new_temp(hh,last_non_nan);
end

% for (hh=1:length(new_temp(:,1)))
%     new_temp(hh,:) = interp1(E_min_MC, new_temp(hh,:), E_min_MC, 'pchip');
% end

[Y,I] = sort(new_temp(:,1), 'ascend');
F = new_temp(I,:);