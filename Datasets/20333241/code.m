%% Chavel~Cullen, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
chavel_cullen_2010.pmid = 20333241;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(chavel_cullen_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/journal.pgen.1000883.s011.xlsx', 'Complete Screen');

% Get the list of ORFs 
hit_orfs = data.raw(:,2);

% And the corresponding data
hit_data = data.raw(:,6);
false_positives = data.raw(:,7);

% Eliminate all white spaces & capitalize
hit_orfs = clean_orf(hit_orfs);

% Eliminate anything that is empty/numeric
inds = find(cellfun(@isnumeric, hit_orfs));
hit_orfs(inds) = [];
hit_data(inds, :) = [];
false_positives(inds, :) = [];

% If in gene name form, transform into ORF name
[hit_orfs, translated, ambiguous] = translate(hit_orfs);

% A couple of manual fixes
hit_orfs(ismember(hit_orfs, {'YLR287-A'})) = {'YLR287C-A'};

inds = find(~is_orf(hit_orfs));
hit_orfs(inds) = [];
hit_data(inds,:) = [];
false_positives(inds,:) = [];

% Eliminate false positives
inds = find(~cellfun(@isnumeric, false_positives));
hit_data(inds) = {NaN};

% Make sure all the data are numbers
inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {NaN};
hit_data = cell2mat(hit_data);

% Average data for identical ORFs that appear multiple times
[hit_orfs,hit_data] = grpstats(hit_data, hit_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
chavel_cullen_2010.orfs = hit_orfs;
chavel_cullen_2010.ph = hit_data_names;
chavel_cullen_2010.data = hit_data;
chavel_cullen_2010.dataset_ids = hit_data_ids;

%% Save

save('./chavel_cullen_2010.mat','chavel_cullen_2010');

%% Print out

fid = fopen('./chavel_cullen_2010.txt','w');
write_matrix_file(fid, chavel_cullen_2010.orfs, chavel_cullen_2010.ph, chavel_cullen_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(chavel_cullen_2010)
end

end
