%% Kanki~Klionsky, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
kanki_klionsky_2009.pmid = 19793921;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(kanki_klionsky_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/TableS2-2.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data
hits_orfs = data.raw(4:end,1);
hits_scores = data.raw(4:end,2);

% Make sure all the indices are numeric
inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds) = [];

% Eliminate all white spaces & capitalize
hits_orfs = clean_orf(hits_orfs);

% If in gene name form, transform into ORF name
[hits_orfs, translated, ambiguous] = translate(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
hits_orfs(inds) = [];
hits_scores(inds) = [];

% Make the data
hits_scores(strcmp('x', hits_scores)) = {NaN};
hits_scores(strcmp('No', hits_scores)) = {-1};
hits_scores(strcmp('Yes', hits_scores)) = {1};
hits_scores = cell2mat(hits_scores);

% If the same strain is present more than once, average its values
[hits_orfs, hits_scores] = grpstats(hits_scores, hits_orfs, {'gname','mean'});
hits_scores(hits_scores == 0) = NaN;

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [160];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
kanki_klionsky_2009.orfs = hits_orfs;
kanki_klionsky_2009.ph = hit_data_names;
kanki_klionsky_2009.data = hits_scores;
kanki_klionsky_2009.dataset_ids = hit_data_ids;

%% Save

save('./kanki_klionsky_2009.mat','kanki_klionsky_2009');

%% Print out

fid = fopen('./kanki_klionsky_2009.txt','w');
write_matrix_file(fid, kanki_klionsky_2009.orfs, kanki_klionsky_2009.ph, kanki_klionsky_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(kanki_klionsky_2009)
end

end
