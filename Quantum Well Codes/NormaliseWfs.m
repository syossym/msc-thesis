function normalised_wfs = NormaliseWfs(z,wfs)

[x,y] = size(wfs);

if (x>y)
   wfs = wfs.'; 
end
normalised_wfs = zeros(size(wfs));
for (ii=1:length(wfs(:,1)))
    dem = sqrt(trapz(z,abs(wfs(ii,:)).^2));
    if (dem~=0)
        normalised_wfs(ii,:) = wfs(ii,:)/sqrt(trapz(z,abs(wfs(ii,:)).^2));
    else
        normalised_wfs(ii,:) = wfs(ii,:);
    end
end
