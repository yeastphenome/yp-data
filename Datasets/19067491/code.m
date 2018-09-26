%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
tyedmers_lindquist_2008.pmid = 19067491;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(tyedmers_lindquist_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/hits.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);

% Get the data itself
hit_data = data(:,2);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% Correct one obvious typo
hit_strains(strcmp(hit_strains, 'RRF')) = {'RRF1'};

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16033];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
tyedmers_lindquist_2008.orfs = hit_strains;
tyedmers_lindquist_2008.ph = hit_data_names;
tyedmers_lindquist_2008.data = hit_data;
tyedmers_lindquist_2008.dataset_ids = hit_data_ids;

%% Save

save('./tyedmers_lindquist_2008.mat','tyedmers_lindquist_2008');

%% Print out

fid = fopen('./tyedmers_lindquist_2008.txt','w');
write_matrix_file(fid, tyedmers_lindquist_2008.orfs, tyedmers_lindquist_2008.ph, tyedmers_lindquist_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(tyedmers_lindquist_2008)
end

end

