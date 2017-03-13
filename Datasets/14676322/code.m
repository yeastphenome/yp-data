%% Warringer~Blomberg, 2003
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

warringer_blomberg_2003.pmid = 14676322;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(warringer_blomberg_2003.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%%

phenotypes = {'Growth, lag phase';'Growth, exponential growth rate';'Growth, saturation level'};
treatments = {'NaCl, 0.85 M'};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/ORIG130305_LPI NaCl.xlsx', 'LPI');

hit_strains = data.raw(5:end,1);
hit_data = cell2mat(data.raw(5:end,2:4));

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [50 49 51]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
warringer_blomberg_2003.orfs = hit_strains;
warringer_blomberg_2003.ph = hit_data_names;
warringer_blomberg_2003.data = hit_data;
warringer_blomberg_2003.dataset_ids = hit_data_ids;

%% Save

save('./warringer_blomberg_2003.mat','warringer_blomberg_2003');

%% Print out

fid = fopen('./warringer_blomberg_2003.txt','w');
write_matrix_file(fid, warringer_blomberg_2003.orfs, warringer_blomberg_2003.ph, warringer_blomberg_2003.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(warringer_blomberg_2003)
end

end

