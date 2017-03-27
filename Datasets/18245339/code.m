%% Abe~Minegishi, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
abe_minegishi_2008.pmid = 18245339;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(abe_minegishi_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested strains

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/mat_alpha_041902.xlsx', 'mat_alpha_041902.txt');

tested_orfs = tested.raw(4:end,3);

tested_orfs = clean_orf(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));  

tested_orfs(inds) = [];

tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Table1_Abe_Genetics.xlsx', 'Table 1 (2)');

hits_orfs = data.raw(7:end,3);
hits_data = data.raw(7:end,[26 31]);

hits_orfs = clean_orf(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));  

hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_data = cell2mat(hits_data);
hits_data = hits_data/100;  % transform percent into fractions

[missing, ix] = setdiff(hits_orfs, tested_orfs);    % nothing missing

% Average data for identical ORFs that appear multiple times
[hits_orfs, hits_data] = grpstats(hits_data, hits_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [537; 538];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
abe_minegishi_2008.orfs = tested_orfs;
abe_minegishi_2008.ph = hit_data_names;
abe_minegishi_2008.data = zeros(length(abe_minegishi_2008.orfs),length(abe_minegishi_2008.ph));
abe_minegishi_2008.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, abe_minegishi_2008.orfs);
abe_minegishi_2008.data(ind2,:) = hits_data(ind1,:);

%% Save

save('./abe_minegishi_2008.mat','abe_minegishi_2008');

%% Print out

fid = fopen('./abe_minegishi_2008.txt','w');
write_matrix_file(fid, abe_minegishi_2008.orfs, abe_minegishi_2008.ph, abe_minegishi_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(abe_minegishi_2008)
end

end
