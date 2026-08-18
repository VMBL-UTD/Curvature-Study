% Calculate effective von Mises stress
% Input: stress_path --> .csv file of stress values output from FEBio sim
% Output: Effective stress values per element for each time step

function [stress] = calculateEffectiveStress(stress_path)

if isnumeric(stress_path)
    data = stress_path;
else
    data = readmatrix(stress_path); 
end

data(isnan(data(:,1)),:) = []; % Removes text
data = sortrows(data,1); % Delete

% Extract stress tensor elements
elem_id = data(:,1);
sx = data(:,2);
sy = data(:,3);
sz = data(:,4);
sxy = data(:,5);
syz = data(:,6);
sxz = data(:,7);

% Effective stress equation (via FEBio Documentation)
stress =  sqrt(0.5.* ((sx - sy).^2 + (sy - sz).^2 + (sx - sz).^2 + 6.*(sxy.^2 + sxz.^2 + syz.^2)));

M = unique(elem_id); % Amount of unique values 
num_steps = abs(size(elem_id,1)/size(M,1)); % Consolidates repetitive numbers into their time steps
stress = reshape(stress,num_steps,[])';

end