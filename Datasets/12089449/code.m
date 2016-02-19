%% Jorgensen~Tyers, 2002
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
jorgensen_tyers_2002.pmid = 12089449;

phenotypes = {'cell size'};
treatments = {''};


%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/data.xlsx', 'haps');

% Get the list of ORFs and the correponding data 
hit_strains = data(:,1);

% Get the data itself
hit_data = data(:,3); 
    
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Normalize data to WT
hit_data(~cellfun(@isnumeric, hit_data)) = {NaN};
hit_data = cell2mat(hit_data);
ind_wt = find(strcmp('WT', hit_strains));
hit_data = hit_data ./ repmat(hit_data(ind_wt,:),length(hit_strains),1);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

%% Prepare final dataset

jorgensen_tyers_2002.orfs = hit_strains;
jorgensen_tyers_2002.ph = strcat(phenotypes, '; ', treatments);
jorgensen_tyers_2002.data = hit_data;

%% Save

save('./jorgensen_tyers_2002.mat','jorgensen_tyers_2002');

end

