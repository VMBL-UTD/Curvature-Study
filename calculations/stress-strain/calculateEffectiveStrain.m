% Calculate effective Lagrange strain
% Input: strain_path --> .csv file of strain values output from FEBio sim
% Output: Effective strain values per element for each time step

function [strain] = calculateEffectiveStrain(strain_path)

if isnumeric(strain_path)
    data = strain_path;
else
    data = readmatrix(strain_path); 
end

data(isnan(data(:,1)),:) = []; % Removes text 
data = sortrows(data,1); % Delete

% Extract strain tensor elements
elem_id = data(:,1);
E11 = data(:,2);
E22 = data(:,3);
E33 = data(:,4);
E12 = data(:,5);
E23 = data(:,6);
E13 = data(:,7);

% Effective strain equation (via FEBio Documentation)
strain = sqrt(E11.^2 + E22.^2 + E33.^2 - E12.*E13 - E12.*E23 - E23.*E13 + 3*(E12.^2 + E13.^2 + E23.^2));

C = unique(elem_id); % Amount of unique values 
num_steps = abs(size(elem_id,1)/size(C,1)); % Consolidates repetitive numbers into their time steps 
strain = reshape(strain,num_steps,[])';

end