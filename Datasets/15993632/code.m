%% Corbacho~Hernandez, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
corbacho_hernandez_2005.pmid = 15993632;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(corbacho_hernandez_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};


%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('textread','./raw_data/Corbacho_tables.txt','%s %s %s %s %*[^\n]','delimiter','\t');

% Get the list of ORFs and the correponding data 
hit_strains = data{1};

% Get the data itself
hit_data = data{4};
    
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

% One phenotype observed only in Mat-alpha strain. Remove the annotation
inds = find(strcmp('1 (alpha)', hit_data));
hit_data(inds) = {'1'};
hit_data = cellfun(@str2num, hit_data);

% Transform data from: 1 = large decrease of staining, 3 = small decrease of staining
% to: -3 = super-low staining, -1 = low staining
hit_data = hit_data-4;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [582];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

corbacho_hernandez_2005.orfs = hit_strains;
corbacho_hernandez_2005.ph = hit_data_names;
corbacho_hernandez_2005.data = hit_data;
corbacho_hernandez_2005.dataset_ids = hit_data_ids;

%% Save

save('./corbacho_hernandez_2005.mat','corbacho_hernandez_2005');

%% Print out

fid = fopen('./corbacho_hernandez_2005.txt','w');
write_matrix_file(fid, corbacho_hernandez_2005.orfs, corbacho_hernandez_2005.ph, corbacho_hernandez_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(corbacho_hernandez_2005)
end

end

