%% Matetic~Smith, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
matecic_smith_2010.pmid = 20421943;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(matecic_smith_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data- downtags

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/journal.pgen.1000921.s002.xlsx', 'DN-tag signal ratio, log2 ranks');

% Get the list of ORFs and the correponding data 
hit_strains_dn = data(5:end,1);

% Get the data itself
hit_data_dn = data(5:end,[11:13 15:17]); 
hit_data_dn = cell2mat(hit_data_dn);

% Eliminate all white spaces & capitalize
hit_strains_dn = clean_orf(hit_strains_dn);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_dn));
hit_strains_dn(inds) = [];
hit_data_dn(inds, :) = [];

%% Load the data- uptags

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/journal.pgen.1000921.s002.xlsx', 'UP-tag signal ratio, log2 ranks');

% Get the list of ORFs and the correponding data 
hit_strains_up = data(5:end,1);

% Get the data itself
hit_data_up = data(5:end,[11:13 15:17]); 
hit_data_up = cell2mat(hit_data_up);

% Eliminate all white spaces & capitalize
hit_strains_up = clean_orf(hit_strains_up);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_up));
hit_strains_up(inds) = [];
hit_data_up(inds, :) = [];

%% Average out up and down tags
hit_strains = [hit_strains_dn; hit_strains_up];
hit_data = [hit_data_dn; hit_data_up];
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [4712; 5354; 5355; 5356; 5357; 5358];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
matecic_smith_2010.orfs = hit_strains;
matecic_smith_2010.ph = hit_data_names;
matecic_smith_2010.data = hit_data;
matecic_smith_2010.dataset_ids = hit_data_ids;

%% Save

save('./matecic_smith_2010.mat','matecic_smith_2010');

%% Print out

fid = fopen('./matecic_smith_2010.txt','w');
write_matrix_file(fid, matecic_smith_2010.orfs, matecic_smith_2010.ph, matecic_smith_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(matecic_smith_2010)
end

end

