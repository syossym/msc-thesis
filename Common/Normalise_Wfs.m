function normWfs = Normalise_Wfs(z,wfs,type)

[x,y] = size(wfs);

if (x>y)
   wfs = wfs.'; 
end
normWfs = zeros(size(wfs));
for (ii=1:length(wfs(:,1)))
    if (strcmp(type, 'Rel'))
        dem = sqrt(trapz(z,abs(wfs(ii,:)).^2));
    elseif (strcmp(type, 'Abs'))
        dem = sqrt(trapz(1:length(z),abs(wfs(ii,:)).^2));
    end
    if (dem~=0)
        normWfs(ii,:) = wfs(ii,:)/dem;
    else
        normWfs(ii,:) = wfs(ii,:);
    end
end
