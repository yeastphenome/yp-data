%% Auesukaree~Harashima, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
auesukaree_harashima_2009.pmid = 19638689;

phenotypes = {'growth'};
treatments = {'EtOH';'MeOH';'propanol';'NaCl';'H2O2';'37C'};

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat alpha_KOset list.xlsx');
tested_orfs = tested.raw(4:end,2);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(upper(clean_orf(tested_orfs)));

% Load data
[FILENAMES{end+1}, data_hits{1}] = read_data('textread','./raw_data/ethanol_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{2}] = read_data('textread','./raw_data/methanol_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{3}] = read_data('textread','./raw_data/propanol_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{4}] = read_data('textread','./raw_data/nacl_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{5}] = read_data('textread','./raw_data/h2o2_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{6}] = read_data('textread','./raw_data/heat_sensitivity_hits.txt', '%s');

data_hits_orfs = cell(size(data_hits));

data_hits_orfs{1} = unique(translate(data_hits{1}));
[missing, ix] = setdiff(data_hits_orfs{1}, tested_orfs);
data_hits_orfs{1}(strcmp('YHR039C-A', data_hits_orfs{1})) = {'YHR039C-B'};

data_hits_orfs{2} = unique(translate(data_hits{2}));
[missing, ix] = setdiff(data_hits_orfs{2}, tested_orfs);

data_hits_orfs{3} = unique(translate(data_hits{3}));
[missing, ix] = setdiff(data_hits_orfs{3}, tested_orfs);

data_hits_orfs{4} = unique(translate(data_hits{4}));
[missing, ix] = setdiff(data_hits_orfs{4}, tested_orfs);

data_hits_orfs{5} = unique(translate(data_hits{5}));
[missing, ix] = setdiff(data_hits_orfs{5}, tested_orfs);


data_hits_orfs{6} = unique(translate(data_hits{6}));
[missing, ix] = setdiff(data_hits_orfs{6}, tested_orfs);
tested_orfs = [tested_orfs; missing];

auesukaree_harashima_2009.orfs = tested_orfs;
auesukaree_harashima_2009.data = zeros(length(tested_orfs),length(treatments));

for i = 1 : length(treatments)
[~,ind1,ind2] = intersect(data_hits_orfs{i}, tested_orfs);
auesukaree_harashima_2009.data(ind2,i) = -1;
end

auesukaree_harashima_2009.ph = strcat(phenotypes, '; ', treatments);

save('./auesukaree_harashima_2009.mat','auesukaree_harashima_2009');

fid = fopen('./auesukaree_harashima_2009.txt','w');
write_matrix_file(fid, auesukaree_harashima_2009.orfs, auesukaree_harashima_2009.ph, auesukaree_harashima_2009.data);
fclose(fid);

end
