%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
anand_payne_2009.pmid = 19116312;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(anand_payne_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/hits.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);

% Get the data itself
hit_data = data(:,3);

tmp = regexp(hit_strains, '?', 'split');
tmp2 = cellfun(@(x) x{1}, tmp, 'UniformOutput', 0);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(tmp2);

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
hit_data_ids = [122];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
anand_payne_2009.orfs = hit_strains;
anand_payne_2009.ph = hit_data_names;
anand_payne_2009.data = hit_data;
anand_payne_2009.dataset_ids = hit_data_ids;

%% Save

save('./anand_payne_2009.mat','anand_payne_2009');

%% Print out

fid = fopen('./anand_payne_2009.txt','w');
write_matrix_file(fid, anand_payne_2009.orfs, anand_payne_2009.ph, anand_payne_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(anand_payne_2009)
end

end

