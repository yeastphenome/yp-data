%% Irwin~Sandmeyer, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
irwin_sandmeyer_2005.pmid = 15837808;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(irwin_sandmeyer_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table1.xlsx');

% Get the list of ORFs and the correponding data 
hit_strains_1 = data(2:end,1);
hit_strains_2 = data(2:65, 2);

% Get the data itself
hit_data_1 = zeros(length(hit_strains_1),1)-1;
hit_data_2 = zeros(length(hit_strains_2),1)+1;

hit_strains = [hit_strains_1; hit_strains_2];
hit_data = [hit_data_1; hit_data_2];
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [182];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
irwin_sandmeyer_2005.orfs = hit_strains;
irwin_sandmeyer_2005.ph = hit_data_names;
irwin_sandmeyer_2005.data = hit_data;
irwin_sandmeyer_2005.dataset_ids = hit_data_ids;

%% Save

save('./irwin_sandmeyer_2005.mat','irwin_sandmeyer_2005');

%% Print out

fid = fopen('./irwin_sandmeyer_2005.txt','w');
write_matrix_file(fid, irwin_sandmeyer_2005.orfs, irwin_sandmeyer_2005.ph, irwin_sandmeyer_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(irwin_sandmeyer_2005)
end

end

