%% Teixeira~Sa-Correia, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
teixeira_sa_correia_2009.pmid = 19633105;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(teixeira_sa_correia_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

tested_orfs = clean_orf(tested_orfs);
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));

tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/TableS1_suplementary_material.xlsx');

hits_genenames = data.raw(7:end,1);
hits_genenames = clean_genename(hits_genenames);
hits_genenames(strcmp('YGL235w', hits_genenames)) = {'YGL235W'};

inds = find(~is_genename(hits_genenames) & ~is_orf(hits_genenames));
hits_genenames(inds) = [];

[hits_orfs, translated] = translate(hits_genenames);
hits_orfs(~translated) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);
hits_orfs(ix) = [];      % 26 ORFs eliminated from the hit list (the list contains hits from both the hap and the het collection, but the list of tested strains is only available for the hap collection)

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [155];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

teixeira_sa_correia_2009.orfs = tested_orfs;
teixeira_sa_correia_2009.data = zeros(length(tested_orfs),1);

[~,~,ind2] = intersect(hits_orfs, tested_orfs);
teixeira_sa_correia_2009.data(ind2) = -1;

teixeira_sa_correia_2009.ph = hit_data_names;
teixeira_sa_correia_2009.dataset_ids = hit_data_ids;

%% Save

save('./teixeira_sa_correia_2009.mat','teixeira_sa_correia_2009');

%% Print out

fid = fopen('./teixeira_sa_correia_2009.txt','w');
write_matrix_file(fid, teixeira_sa_correia_2009.orfs, teixeira_sa_correia_2009.ph, teixeira_sa_correia_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(teixeira_sa_correia_2009)
end

end
