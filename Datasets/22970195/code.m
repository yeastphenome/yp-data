%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ayer_perrone_2012.pmid = 22970195;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ayer_perrone_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/journal.pone.0044278.s004.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);
hit_strains{1} = 'WT';

% Get the data itself
hit_data = data(2:end,4);

% Eliminate the NaNs
inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

r = '([0-9]{3})(\s)*±([0-9]{1})';
hit_data2 = nan(size(hit_data));
for i = 1 : length(hit_data)
    [tokens, ~] = regexp(hit_data{i}, r, 'tokens','match');
    hit_data2(i) = str2num(tokens{1}{1});
end

% Normalize by WT
hit_data2 = hit_data2 ./ hit_data2(1);
hit_data2(1) = [];
hit_strains(1) = [];

hit_data = hit_data2;

  
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [514];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
ayer_perrone_2012.orfs = hit_strains;
ayer_perrone_2012.ph = hit_data_names;
ayer_perrone_2012.data = hit_data;
ayer_perrone_2012.dataset_ids = hit_data_ids;

%% Save

save('./ayer_perrone_2012.mat','ayer_perrone_2012');

%% Print out

fid = fopen('./ayer_perrone_2012.txt','w');
write_matrix_file(fid, ayer_perrone_2012.orfs, ayer_perrone_2012.ph, ayer_perrone_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ayer_perrone_2012)
end

end

