%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
witkin_cohen_fix_2012.pmid = 22658600;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(witkin_cohen_fix_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS1.xlsx', 'Table 1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,2);

% Get the data itself
hit_data = ones(size(hit_strains));

inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11848];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
witkin_cohen_fix_2012.orfs = hit_strains;
witkin_cohen_fix_2012.ph = hit_data_names;
witkin_cohen_fix_2012.data = hit_data;
witkin_cohen_fix_2012.dataset_ids = hit_data_ids;

%% Save

save('./witkin_cohen_fix_2012.mat','witkin_cohen_fix_2012');

%% Print out

fid = fopen('./witkin_cohen_fix_2012.txt','w');
write_matrix_file(fid, witkin_cohen_fix_2012.orfs, witkin_cohen_fix_2012.ph, witkin_cohen_fix_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(witkin_cohen_fix_2012)
end

end
