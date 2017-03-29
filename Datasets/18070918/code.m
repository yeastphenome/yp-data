%% Ju~Xie, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ju_xie_2008.pmid = 18070918;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ju_xie_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/mat_a_041902.xlsx', 'mat_a_041902');
tested_orfs = data.raw(3:end,2);

tested_orfs = clean_orf(tested_orfs);

tested_orfs(strcmp('YLR287-A', tested_orfs)) = {'YLR287C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));  

tested_orfs(inds) = [];

tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES, genenames] = read_data('textscan', './raw_data/hits.txt', '%s');

orfs = translate(genenames);
data = ones(size(orfs));

missing_orfs = setdiff(orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [166];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
ju_xie_2008.orfs = tested_orfs;
ju_xie_2008.ph = hit_data_names;
ju_xie_2008.data = zeros(length(ju_xie_2008.orfs),length(ju_xie_2008.ph));
ju_xie_2008.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(orfs, ju_xie_2008.orfs);
ju_xie_2008.data(ind2,:) = data(ind1,:);

%% Save

save('./ju_xie_2008.mat','ju_xie_2008');

%% Print out

fid = fopen('./ju_xie_2008.txt','w');
write_matrix_file(fid, ju_xie_2008.orfs, ju_xie_2008.ph, ju_xie_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ju_xie_2008)
end

end
