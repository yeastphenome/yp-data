%% Teixeira~Sa-Correia, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
teixeira_sa_correia_2010.pmid = 20210661;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(teixeira_sa_correia_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/List of strains tested.xlsx');

% Get the list of ORFs
tested_orfs = tested.raw(2:end,1);

% Eliminate all white spaces & capitalize
tested_orfs = clean_orf(tested_orfs);

% If in gene name form, transform into ORF name
[tested_orfs, translated, ambiguous] = translate(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

% Finally, take the unique set
tested_orfs = unique(tested_orfs);

%% Load hit data

[FILENAMES{end+1}, hits_genenames] = read_data('textread','./raw_data/hits_genenames.txt', '%s');

% Eliminate all white spaces & capitalize
hits_genenames = clean_genename(hits_genenames);

% If in gene name form, transform into ORF name
[hits_orfs, translated, ambiguous] = translate(hits_genenames);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
hits_orfs(inds) = [];

% Make sure the that all the hits are part of the tested set
[missing, ix] = setdiff(hits_orfs, tested_orfs);
tested_orfs = [tested_orfs; missing];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [151];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
teixeira_sa_correia_2010.orfs = tested_orfs;
teixeira_sa_correia_2010.ph = hit_data_names;
teixeira_sa_correia_2010.data = zeros(length(teixeira_sa_correia_2010.orfs),length(teixeira_sa_correia_2010.ph));
teixeira_sa_correia_2010.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, teixeira_sa_correia_2010.orfs);
teixeira_sa_correia_2010.data(ind2) = -1;

%% Save

save('./teixeira_sa_correia_2010.mat','teixeira_sa_correia_2010');

%% Print out

fid = fopen('./teixeira_sa_correia_2010.txt','w');
write_matrix_file(fid, teixeira_sa_correia_2010.orfs, teixeira_sa_correia_2010.ph, teixeira_sa_correia_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(teixeira_sa_correia_2010)
end

end
