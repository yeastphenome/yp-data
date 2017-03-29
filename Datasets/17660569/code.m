%% Wilson~van Hoof, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
wilson_vanhoof_2007.pmid = 17660569;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(wilson_vanhoof_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/hit_list.xlsx', 'Sheet1');

% Get the list of ORFs
hit_strains = data(3:end,1);

% Get the data itself
hit_data = data(3:end,3);
    
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds) = [];

hit_data = cellfun(@length, hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [169];

%% Load tested strains

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/deletion database.xlsx', 'Sheet1');

% Get the list of ORFs
tested_strains = data(2:end,2);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% Fix typos
tested_strains(strcmp('YPL006', tested_strains)) = {'YPL006W'};
tested_strains(strcmp('YDR007', tested_strains)) = {'YDR007W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds)); 

tested_strains(inds) = [];

tested_strains = unique(tested_strains);

[missing,~] = setdiff(hit_strains, tested_strains); % 0 found.

%% Prepare final dataset

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
wilson_vanhoof_2007.orfs = tested_strains;
wilson_vanhoof_2007.ph = hit_data_names;
wilson_vanhoof_2007.data = zeros(length(wilson_vanhoof_2007.orfs),length(wilson_vanhoof_2007.ph));
wilson_vanhoof_2007.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, wilson_vanhoof_2007.orfs);
wilson_vanhoof_2007.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./wilson_vanhoof_2007.mat','wilson_vanhoof_2007');

%% Print out

fid = fopen('./wilson_vanhoof_2007.txt','w');
write_matrix_file(fid, wilson_vanhoof_2007.orfs, wilson_vanhoof_2007.ph, wilson_vanhoof_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(wilson_vanhoof_2007)
end

end

