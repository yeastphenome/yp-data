%% Kemmer~Roberge, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
kemmer_roberge_2009.pmid = 19144191;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(kemmer_roberge_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/haploid set.xlsx', 'haploid set');
tested_orfs = tested.raw(6:end,2);

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = clean_orf(tested_orfs);
tested_orfs(strcmp('YYKL138C', tested_orfs)) = {'YKL138C'};

tested_orfs = translate(tested_orfs);

tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, hits_genenames] = read_data('textread','./raw_data/hits_genenames.txt', '%s');

hits_genenames = clean_genename(hits_genenames);
hits_orfs = translate(hits_genenames);

inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));

hits_orfs = unique(hits_orfs);

hits_scores = -ones(length(hits_orfs),1);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [159];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
kemmer_roberge_2009.orfs = tested_orfs;
kemmer_roberge_2009.ph = hit_data_names;
kemmer_roberge_2009.data = zeros(length(kemmer_roberge_2009.orfs),length(kemmer_roberge_2009.ph));
kemmer_roberge_2009.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, kemmer_roberge_2009.orfs);
kemmer_roberge_2009.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./kemmer_roberge_2009.mat','kemmer_roberge_2009');

%% Print out

fid = fopen('./kemmer_roberge_2009.txt','w');
write_matrix_file(fid, kemmer_roberge_2009.orfs, kemmer_roberge_2009.ph, kemmer_roberge_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(kemmer_roberge_2009)
end

end
