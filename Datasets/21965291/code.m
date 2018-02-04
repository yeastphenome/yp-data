%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
zakrzewska_smits_2011.pmid = 21965291;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(zakrzewska_smits_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data 1

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/mc-E10-08-0721-s06.xlsx', 'growth rates');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data1(2:end,1);

% Get the data itself
hit_data1 = data1(2:end,2:3); % if the dataset is discrete or binary
   
% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

inds = find(~cellfun(@isnumeric, hit_data1));
hit_data1(inds) = {NaN};

hit_data1 = cell2mat(hit_data1);

% If the same strain is present more than once, average its values
[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids1 = [16128,16129]';

%% Load the data 2

[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/mc-E10-08-0721-s06.xlsx', 'survival % 95% CI');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains2 = data2(2:end,1);

% Get the data itself
hit_data2 = data2(2:end,2:4:end); % if the dataset is discrete or binary
   
% Eliminate all white spaces & capitalize
hit_strains2 = clean_orf(hit_strains2);

% If in gene name form, transform into ORF name
hit_strains2 = translate(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

hit_data2 = cell2mat(hit_data2);

% If the same strain is present more than once, average its values
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids2 = [16130:16135]';

%% Merge the 2 datasets

hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = nan(length(hit_strains), 8);

[~,ind1,ind2] = intersect(hit_strains, hit_strains1);
hit_data(ind1,1:2) = hit_data1(ind2,:);

[~,ind1,ind2] = intersect(hit_strains, hit_strains2);
hit_data(ind1,3:end) = hit_data2(ind2,:);

hit_data_ids = [hit_data_ids1; hit_data_ids2];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
zakrzewska_smits_2011.orfs = hit_strains;
zakrzewska_smits_2011.ph = hit_data_names;
zakrzewska_smits_2011.data = hit_data;
zakrzewska_smits_2011.dataset_ids = hit_data_ids;

%% Save

save('./zakrzewska_smits_2011.mat','zakrzewska_smits_2011');

%% Print out

fid = fopen('./zakrzewska_smits_2011.txt','w');
write_matrix_file(fid, zakrzewska_smits_2011.orfs, zakrzewska_smits_2011.ph, zakrzewska_smits_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(zakrzewska_smits_2011)
end

end

