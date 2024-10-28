function [MinT,MinDist]=linemindist_t(A1,B1,D1,E1,T1,A2,B2,D2,E2,T2);
%------------------------------------------------------------------------
% linemindist_t function        Given two 2D linear trajectories
%                             calculate the time in which the distance
%                             between the two trajectories is minimal 
%                             and the distance at that time.
%                             Each line is defined by:
%                             x_i = A_i + B_i*(t - T_i)
%                             y_i = D_i + E_i*(t - T_t)
% Input  : - A1
%          - B1
%          - D1
%          - E1
%          - T1
%          - A2
%          - B2
%          - D2
%          - E2
%          - T2
% Output : - Time in which the distance is minimum.
%          - Minimum distance.
% Tested : Matlab 5.3
%     By : Eran O. Ofek         June 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%------------------------------------------------------------------------

MinT = (-A1.*B1+A1.*B2+B1.^2.*T1-B1.*T1.*B2+A2.*B1-A2.*B2-B2.*T2.*B1+B2.^2.*T2-D1.*E1+D1.*E2+E1.^2.*T1-E1.*T1.*E2+D2.*E1-D2.*E2-E2.*T2.*E1+E2.^2.*T2)./(B1.^2-2.*B1.*B2+B2.^2+E1.^2-2.*E1.*E2+E2.^2);

%--- this come from ---
%diff('(a1+b1*(t-T1) - a2-b2*(t-T2))^2+(d1+e1*(t-T1) - d2-e2*(t-T2))^2','t')
%vectorize(solve('2*(a1+b1*(t-T1)-a2-b2*(t-T2))*(b1-b2)+2*(d1+e1*(t-T1)-d2-e2*(t-T2))*(e1-e2)=0','t'))

MinDist = sqrt((A1 + B1.*(MinT - T1) - A2 - B2.*(MinT - T2)).^2 + (D1 + E1.*(MinT - T1) - D2 - E2.*(MinT - T2)).^2);
