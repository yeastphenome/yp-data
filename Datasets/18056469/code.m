%% Xia~Flores-Rozas, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
xia_flores_rozas_2007.pmid = 18056469;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(xia_flores_rozas_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat_a.xlsx', 'mat_a_041902');
tested_orfs = tested.raw(3:end,2);

tested_orfs = clean_orf(tested_orfs);

tested_orfs(strcmp('YLR287-A', tested_orfs)) = {'YLR287C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));  

tested_orfs(inds) = [];
tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, DATA] = read_data('textread','./raw_data/hits_genenames.txt', '%s %s', 'delimiter', '\t');

hits_genenames = DATA{1};
hits_scores_txt = DATA{2};

hits_genenames = clean_genename(hits_genenames);
hits_orfs = translate(hits_genenames);
hits_scores = -cellfun(@length, hits_scores_txt);

hits_orfs(strcmp('TFP3', hits_orfs)) = {'YPL234C'};

[hits_orfs, hits_scores] = grpstats(hits_scores, hits_orfs, {'gname','mean'});

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [97];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
xia_flores_rozas_2007.orfs = tested_orfs;
xia_flores_rozas_2007.ph = hit_data_names;
xia_flores_rozas_2007.data = zeros(length(xia_flores_rozas_2007.orfs),length(xia_flores_rozas_2007.ph));
xia_flores_rozas_2007.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, xia_flores_rozas_2007.orfs);
xia_flores_rozas_2007.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./xia_flores_rozas_2007.mat','xia_flores_rozas_2007');

%% Print out

fid = fopen('./xia_flores_rozas_2007.txt','w');
write_matrix_file(fid, xia_flores_rozas_2007.orfs, xia_flores_rozas_2007.ph, xia_flores_rozas_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(xia_flores_rozas_2007)
end

end
