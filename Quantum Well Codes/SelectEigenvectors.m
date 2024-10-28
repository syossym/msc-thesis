function [D_f, V_f] = SelectEigenvectors(D,V)

index = [];

for (ii=1:2:length(V(1,:))-1)
    if (sum(abs(V(:,ii))) == 0 && sum(abs(V(:,ii+1))) ~= 0)
        index = [index, ii+1];
        continue;
    elseif (sum(abs(V(:,ii))) ~= 0 && sum(abs(V(:,ii+1))) == 0)
        index = [index, ii];
        continue;
    elseif (sum(abs(V(:,ii))) == 0 && sum(abs(V(:,ii+1))) == 0)
        continue;
    end

    if(sum(abs(V(:,ii)) > sum(abs(V(:,ii+1)))))
        index = [index, ii];
    else
        index = [index, ii+1];
    end
end

V_f = V(:,index);
D_f = D(index);