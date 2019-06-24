%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
henriques_sa_correia_2017.pmid = 28428821;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(henriques_sa_correia_2017.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/13068_2017_781_MOESM1_ESM.xlsx', 'Table S1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data1(4:end,1);

% Get the data itself
hit_data = data1(4:end,3);

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

hit_data = cellfun(@length, hit_data);
hit_data = -hit_data;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

%% Resistant genes

[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/13068_2017_781_MOESM2_ESM.xlsx', 'Table S2');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains2 = data2(4:end,1);

% Get the data itself
hit_data2 = ones(length(hit_strains2),1);

inds = find(cellfun(@isnumeric, hit_strains2));
hit_strains2(inds) = [];
hit_data2(inds,:) = [];
   
% Eliminate all white spaces & capitalize
hit_strains2 = clean_genename(hit_strains2);

% If in gene name form, transform into ORF name
hit_strains2 = translate(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

hit_strains2(inds) = [];
hit_data2(inds,:) = [];

% If the same strain is present more than once, average its values
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

%% Merge

hit_strains = [hit_strains; hit_strains2];
hit_data = [hit_data; hit_data2];

%%
% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16264];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('xlsread','./raw_data/List of strains tested.xlsx', 'Tabelle2');

tested_strains = tested_strains(:,1);

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

% If it seems reasonable, add the missing hits to the list of tested strains
tested_strains = [tested_strains; missing];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
henriques_sa_correia_2017.orfs = tested_strains;
henriques_sa_correia_2017.ph = hit_data_names;
henriques_sa_correia_2017.data = zeros(length(henriques_sa_correia_2017.orfs),length(henriques_sa_correia_2017.ph));
henriques_sa_correia_2017.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, henriques_sa_correia_2017.orfs);
henriques_sa_correia_2017.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./henriques_sa_correia_2017.mat','henriques_sa_correia_2017');

%% Print out

fid = fopen('./henriques_sa_correia_2017.txt','w');
write_matrix_file(fid, henriques_sa_correia_2017.orfs, henriques_sa_correia_2017.ph, henriques_sa_correia_2017.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(henriques_sa_correia_2017)
end

end

