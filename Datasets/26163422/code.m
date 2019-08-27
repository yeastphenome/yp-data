%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
yimit_brown_2015.pmid = 26163422;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(yimit_brown_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS4.xlsx', 'Table S4');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(4:50,2);

% Get the data itself
hit_data_t = data(4:50,3);
hit_data = zeros(length(hit_data_t),2);
hit_data(strcmp('Phleo treated cells have fewer foci', hit_data_t),1) = -1;
hit_data(strcmp('Untreated cells have more foci', hit_data_t),2) = 1;
   
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
hit_data_ids = [16255; 16320];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('readtable','./raw_data/FG_array_genes.txt', 'ReadVariableNames', 0);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains.Var1);

% If in gene name form, transform into ORF name
tested_strains = translate(tested_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

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
yimit_brown_2015.orfs = tested_strains;
yimit_brown_2015.ph = hit_data_names;
yimit_brown_2015.data = zeros(length(yimit_brown_2015.orfs),length(yimit_brown_2015.ph));
yimit_brown_2015.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, yimit_brown_2015.orfs);
yimit_brown_2015.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./yimit_brown_2015.mat','yimit_brown_2015');

%% Print out

fid = fopen('./yimit_brown_2015.txt','w');
write_matrix_file(fid, yimit_brown_2015.orfs, yimit_brown_2015.ph, yimit_brown_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(yimit_brown_2015)
end

end

