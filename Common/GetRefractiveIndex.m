function n = GetRefractiveIndex(material, E, x)

load Refractive_Index_Fit_Params.mat;

switch (material)
    case 'GaAs',
        %n = polyval(re_n_GaAs_fit, E);
        n = interp1(re_n_GaAs(:,2), re_n_GaAs(:,3), E, 'cubic');
        %n = n + 1i*interp1(im_n_GaAs(:,2), im_n_GaAs(:,3), E, 'cubic')
        E_im = E(E>=1.526);
        if (~isempty(E_im))
            %n = n + 1i*polyval(im_n_GaAs_fit, E);
            im_n = interp1(im_n_GaAs(:,2), im_n_GaAs(:,3), E_im, 'cubic');
            n = n + 1i*padarray(im_n.', length(E)-length(E_im), im_n(1), 'pre').';
        end
    case 'AlAs',
        %n = polyval(re_n_AlAs_fit, E);
        n = interp1(re_n_AlAs(:,2), re_n_AlAs(:,3), E, 'cubic');
        n_im = 1.7e-3;
        n = n + 1i*n_im;
    case 'GaAlAs',
        if (nargin == 2)
            %n = polyval(re_n_GaAlAs_0_1_fit, E);
            n = interp1(re_n_AlGaAs_0_1(:,2), re_n_AlGaAs_0_1(:,3), E, 'cubic');
            %n = n + 1i*interp1(im_n_AlGaAs_0_1(:,2), im_n_AlGaAs_0_1(:,3), E, 'cubic');
            E_im = E(E>=1.637);
            if (~isempty(E_im))
                %n = n + 1i*polyval(im_n_GaAs_fit, E);
                im_n = interp1(im_n_AlGaAs_0_1(:,2), im_n_AlGaAs_0_1(:,3), E_im, 'cubic');
                n = n + 1i*padarray(im_n.', length(E)-length(E_im), im_n(1), 'pre').';
            end
        else
            n_GaAs = GetRefractiveIndex('GaAs', E);
            n_AlAs = GetRefractiveIndex('AlAs', E);
            n = x*n_AlAs + (1-x)*n_GaAs;
        end
end