%% Dakshinamurthy~Garfinkel, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
dakshinamurthy_garfinkel_2010.pmid = 20498295;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(dakshinamurthy_garfinkel_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Matalphakos counted.xlsx', 'Sheet2');

% Get the list of ORFs and the correponding data 
tested_orfs = tested.raw(6:end,2);

% Eliminate all empty spaces, or numeric spaces
inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

% Eliminate all white spaces & capitalize
tested_orfs = clean_orf(tested_orfs);

% Translate
[tested_orfs, translated, ambiguous] = translate(tested_orfs);

% Fix typos
tested_orfs(find(strcmp('YLR228', tested_orfs))) = {'YLR228C'};
tested_orfs(find(strcmp('YMR062', tested_orfs))) = {'YMR062C'};
tested_orfs(find(strcmp('YYKL138C', tested_orfs))) = {'YKL138C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/TableS1-2.xlsx', 'Sheet2');

% Get the list of ORFs and the correponding data 
hits_orfs = data.raw(:,2);

% Get the data itself
hits_data = data.raw(:,3);

% Eliminate all empty spaces, or numeric spaces
inds = find(cellfun(@isempty, hits_orfs) | cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

% Eliminate all white spaces & capitalize
hits_orfs = clean_orf(hits_orfs);

% Translate
[hits_orfs, translated, ambiguous] = translate(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

% Transform dat to numeric
hits_data = cellfun(@length, hits_data);
hits_data = -hits_data; % these data describes a decrease in the number if His+ cells, thus a decrease in the frequency of Ty1 transposition

% Make sure the that all the hits are part of the tested set
missing = setdiff(hits_orfs, tested_orfs);

% Adding 4 ORFs to the list of tested
tested_orfs = [tested_orfs; missing];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [102];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
dakshinamurthy_garfinkel_2010.orfs = tested_orfs;
dakshinamurthy_garfinkel_2010.ph = hit_data_names;
dakshinamurthy_garfinkel_2010.data = zeros(length(dakshinamurthy_garfinkel_2010.orfs),length(dakshinamurthy_garfinkel_2010.ph));
dakshinamurthy_garfinkel_2010.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, dakshinamurthy_garfinkel_2010.orfs);
dakshinamurthy_garfinkel_2010.data(ind2,:) = hits_data(ind1,:);

%% Save

save('./dakshinamurthy_garfinkel_2010.mat','dakshinamurthy_garfinkel_2010');

%% Print out

fid = fopen('./dakshinamurthy_garfinkel_2010.txt','w');
write_matrix_file(fid, dakshinamurthy_garfinkel_2010.orfs, dakshinamurthy_garfinkel_2010.ph, dakshinamurthy_garfinkel_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(dakshinamurthy_garfinkel_2010)
end

end
