%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
elbaz_alon_schuldiner_2014.pmid = 25026036;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(elbaz_alon_schuldiner_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/mmc2.xlsx', 'Table S1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(3:end,1);

% Get the data itself
hit_data = data(3:end,4); % if the dataset is discrete or binary

inds = find(cellfun(@isnumeric, hit_data));
hit_data(inds) = [];
hit_strains(inds) = [];

[phenotypes, ia, ib] = unique(hit_data);
phenotypes_severity = [-1,-1,1,1,-2];

hit_data2 = phenotypes_severity(ib);
hit_data = hit_data2';
   
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
hit_data_ids = [16000];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains_data] = read_data('xlsread','./raw_data/KO_DAmP_ORFs.xlsx', 'Sheet1');

tested_strains = tested_strains_data(3:end,1);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If in gene name form, transform into ORF name
tested_strains = translate(tested_strains);

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'YOLO57W'})) = {'YOL057W'};
tested_strains(ismember(tested_strains, {'YOLO62C'})) = {'YOL062C'};
tested_strains(ismember(tested_strains, {'YBRF182C-A'})) = {'YBR182C-A'};
tested_strains(ismember(tested_strains, {'YLR287-A'})) = {'YLR287C-A'};
tested_strains(ismember(tested_strains, {'YJL206-A'})) = {'YJL206C-A'};

% If not possible, eliminate the entry
tested_strains(ismember(tested_strains, {''})) = [];

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,inds] = setdiff(hit_strains, tested_strains);
disp(missing);

% Removing the missing strains (they were tested as DAMP strains, not
% deletions?)
hit_strains(inds) = [];
hit_data(inds) = [];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
elbaz_alon_schuldiner_2014.orfs = tested_strains;
elbaz_alon_schuldiner_2014.ph = hit_data_names;
elbaz_alon_schuldiner_2014.data = zeros(length(elbaz_alon_schuldiner_2014.orfs),length(elbaz_alon_schuldiner_2014.ph));
elbaz_alon_schuldiner_2014.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, elbaz_alon_schuldiner_2014.orfs);
elbaz_alon_schuldiner_2014.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./elbaz_alon_schuldiner_2014.mat','elbaz_alon_schuldiner_2014');

%% Print out

fid = fopen('./elbaz_alon_schuldiner_2014.txt','w');
write_matrix_file(fid, elbaz_alon_schuldiner_2014.orfs, elbaz_alon_schuldiner_2014.ph, elbaz_alon_schuldiner_2014.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(elbaz_alon_schuldiner_2014)
end

end
