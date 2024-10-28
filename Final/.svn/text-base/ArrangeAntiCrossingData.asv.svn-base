function f = ArrangeAntiCrossingData(E_r_m_MC,E_r_m_MCQW)

% Remove the outliers
for (ii=1:length(E_r_m_MCQW(:,1)))
    f{ii}.E_MC = E_r_m_MC(E_r_m_MCQW(ii,:)~=0);
    f{ii}.E_MCQW = E_r_m_MCQW(ii, E_r_m_MCQW(ii,:)~=0);
    f{ii}.Max = max(f{ii}.E_MCQW);
    f{ii}.Min = min(f{ii}.E_MCQW);
end

