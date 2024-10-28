% Try to run this file to find problem
% To be inserted in the folder of the OSCAR file
% works for 2 mirror cavities (with the default name for the mirrors)

CreateField
Find_error = false;
% test the mirrors profile

if max(max(abs(imag(Mirror.ITM_cav)))) ~= 0
    disp('Incorrect definition for the ITM profile (reflective side)')
    disp('Check the matrix Mirror.ITM_cav')
    Find_error = true;
end

if max(max(abs(imag(Mirror.ETM_cav)))) ~= 0
    disp('Incorrect definition for the ETM profile (reflective side)')
    disp('Check the matrix Mirror.ETM_cav')
    Find_error = true;
end

if max(max(isnan(Field.Start))) ~= 0
    disp('Incorrect definition for the input laser beam')
    disp('Check the matrix Field.Start')
    Find_error = true;
end

if exist('Mirror.ITM_trans','var')
    if max(max(isnan(Mirror.ITM_trans))) ~= 0
        disp('Problem with the matrix Mirror.ITM_trans')
        Find_error = true;
    end
end

if exist('Mirror.ETM_trans','var')
    if max(max(isnan(Mirror.ETM_trans))) ~= 0
        disp('Problem with the matrix Mirror.ETM_trans')
        Find_error = true;
    end
end

if exist('Mirror.ITM_ref','var')
    if max(max(isnan(Mirror.ITM_ref))) ~= 0
        disp('Problem with the matrix Mirror.ITM_ref')
        Find_error = true;
    end
end


if ITM.R > 1
     disp('Reflectivity ITM >1')
     Find_error = true;
end

if ETM.R > 1
     disp('Reflectivity ETM >1')
     Find_error = true;
end

% Check the resolution of the mirror grid
if (max(max(diff(Mirror.ITM_cav))) > Laser.lambda) || (max(max(diff(Mirror.ETM_cav))) > Laser.lambda)
   disp('It may be wise to increase the grid resolution')
   Find_error = true;
end

% Calculate the cavity g factor (no mirror map)
g_cav = (1 - Length_cav/ITM.RofC)*(1 - Length_cav/ETM.RofC);
if (g_cav < 0) && (g_cav > 1)
   disp('Cavity without mirror map is unstable')
   Find_error = true;
end
    
    
    
% Run Cavity stability test
if ~Find_error
    Length.nb_propa_field = 100;
    Propagate_Field

    RT_loss = 1 - (ITM.R*ETM.R);

    for m = 1:(Length.nb_propa_field-1)
        Ploss = 1-Calculate_power(Field.propa(:,:,m+1))/(Calculate_power(Field.propa(:,:,m)));
        if (Ploss /RT_loss) > 2
            disp('Maybe unstable cavity or large diffraction loss')
            Find_error = true;
            break
        end
    end
end

if ~Find_error
    disp('After a quick glance, everything looks fine')
end