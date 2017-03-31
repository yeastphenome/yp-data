%% Mira~Sa-Correia, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
mira_sa_correia_2010.pmid = 20973990;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(mira_sa_correia_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/List of strains tested.xlsx');

% Get the list of ORFs and the correponding data 
tested_orfs = tested.raw(2:end,1);

% Eliminate all white spaces & capitalize
tested_orfs = clean_orf(tested_orfs);

% If in gene name form, transform into ORF name
tested_orfs = translate(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

% Finally, take the unique set
tested_orfs = unique(tested_orfs);

%% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/1475-2859-9-79-s1.xlsx');

% Get the list of ORFs and the correponding data 
hits_genenames = data.raw(9:end,1);
hits_data = data.raw(9:end,3);

% Remove empty gene names
inds = find(cellfun(@isempty, hits_genenames) | cellfun(@isnumeric, hits_genenames));
hits_genenames(inds) = [];
hits_data(inds,:) = [];

% Eliminate all white spaces & capitalize
hits_genenames = clean_genename(hits_genenames);

% If in gene name form, transform into ORF name
[hits_orfs, translated, ambiguous] = translate(hits_genenames);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

% Translate data from symbols to numbers
hits_data(strcmp('++', hits_data)) = {-2};
hits_data(strcmp('+', hits_data)) = {-1};

% Transform into a numeric array
hits_data = cell2mat(hits_data);
hits_data(isnan(hits_data)) = 0;

% Add back the ORFs found in the hits and not in the tested
[missing, ix] = setdiff(hits_orfs, tested_orfs);    % 9 ORFs
tested_orfs = [tested_orfs; missing];

% Average out any repeated data
[hits_data,hits_orfs] = grpstats(hits_data, hits_orfs, {'mean','gname'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [101];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
mira_sa_correia_2010.orfs = tested_orfs;
mira_sa_correia_2010.ph = hit_data_names;
mira_sa_correia_2010.data = zeros(length(mira_sa_correia_2010.orfs),length(mira_sa_correia_2010.ph));
mira_sa_correia_2010.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
mira_sa_correia_2010.data(ind2) = hits_data(ind1);

%% Save
save('./mira_sa_correia_2010.mat','mira_sa_correia_2010');

%% Print out
fid = fopen('./mira_sa_correia_2010.txt','w');
write_matrix_file(fid, mira_sa_correia_2010.orfs, mira_sa_correia_2010.ph, mira_sa_correia_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(mira_sa_correia_2010)
end

end
