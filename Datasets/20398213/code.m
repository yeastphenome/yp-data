%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ruiz_roig_de_nadal_2010.pmid = 20398213;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ruiz_roig_de_nadal_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS1.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);
hit_data = cell2mat(data(:,4));
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

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
hit_data_ids = [15983];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('xlsread','./raw_data/KO_386.xlsx', 'Sheet1');

tested_strains = tested_strains(:,7);

inds = find(cellfun(@isnumeric, tested_strains));
tested_strains(inds) = [];

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If in gene name form, transform into ORF name
tested_strains = translate(tested_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

tested_strains(inds) = [];

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
ruiz_roig_de_nadal_2010.orfs = tested_strains;
ruiz_roig_de_nadal_2010.ph = hit_data_names;
ruiz_roig_de_nadal_2010.data = zeros(length(ruiz_roig_de_nadal_2010.orfs),length(ruiz_roig_de_nadal_2010.ph));
ruiz_roig_de_nadal_2010.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, ruiz_roig_de_nadal_2010.orfs);
ruiz_roig_de_nadal_2010.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./ruiz_roig_de_nadal_2010.mat','ruiz_roig_de_nadal_2010');

%% Print out

fid = fopen('./ruiz_roig_de_nadal_2010.txt','w');
write_matrix_file(fid, ruiz_roig_de_nadal_2010.orfs, ruiz_roig_de_nadal_2010.ph, ruiz_roig_de_nadal_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ruiz_roig_de_nadal_2010)
end

end

