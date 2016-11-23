%% Gonzalez~Valero, 2016
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
gonzalez_valero_2016.pmid = 27733850;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(gonzalez_valero_2016.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the first dataset: sorbitol

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Data_Sheet_1.xlsx', 'Sorbitol vs. Basal');

% Get the list of ORFs and the correponding data 
hit_strains = data(2:end,1);

% Get the data itself
hit_data = data(2:end,2);
hit_data = cell2mat(hit_data);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Fix strains with multiple replications
indx = find(~cellfun(@isempty, strfind(hit_strains,'-3')));
indx1 = find(~cellfun(@isempty, strfind(hit_strains,'-2')));
indx2 = find(~cellfun(@isempty, strfind(hit_strains,'-4')));
indx = [indx; indx1; indx2];
hit_strains(indx) = cellfun(@(x) x(1:end-2), hit_strains(indx), 'un', 0);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YPL183W-A4'})) = {'YPL183W-A'};
hit_strains(ismember(hit_strains, {'YPL183W-A3'})) = {'YPL183W-A'};
hit_strains(ismember(hit_strains, {'YOR298C-A2'})) = {'YOR298C-A'};
hit_strains(ismember(hit_strains, {'YPL183W-A2'})) = {'YPL183W-A'};
hit_strains(ismember(hit_strains, {'YOR298C-A4'})) = {'YOR298C-A'};
hit_strains(ismember(hit_strains, {'YOR298C-A3'})) = {'YOR298C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

%% Second dataset: Glucose

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Data_Sheet_1.xlsx', 'Glucose vs. Basal');

% Get the list of ORFs and the correponding data 
hit_strains2 = data(2:end,1);

% Get the data itself
hit_data2 = data(2:end,2);
hit_data2 = cell2mat(hit_data2);
   
% Eliminate all white spaces & capitalize
hit_strains2 = clean_orf(hit_strains2);

% Fix strains with multiple replications
indx = find(~cellfun(@isempty, strfind(hit_strains2,'-3')));
indx1 = find(~cellfun(@isempty, strfind(hit_strains2,'-2')));
indx2 = find(~cellfun(@isempty, strfind(hit_strains2,'-4')));
indx = [indx; indx1; indx2];
hit_strains2(indx) = cellfun(@(x) x(1:end-2), hit_strains2(indx), 'un', 0);

% If possible, fix the problem (typos, omissions etc.)
hit_strains2(ismember(hit_strains2, {'YPL183W-A4'})) = {'YPL183W-A'};
hit_strains2(ismember(hit_strains, {'YPL183W-A3'})) = {'YPL183W-A'};
hit_strains(ismember(hit_strains, {'YOR298C-A2'})) = {'YOR298C-A'};
hit_strains(ismember(hit_strains, {'YPL183W-A2'})) = {'YPL183W-A'};
hit_strains(ismember(hit_strains, {'YOR298C-A4'})) = {'YOR298C-A'};
hit_strains(ismember(hit_strains, {'YOR298C-A3'})) = {'YOR298C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));

% If the same strain is present more than once, average its values
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

%% Combine data from two datasets

all_strains =[hit_strains; hit_strains2];
all_strains = unique(all_strains);

all_data = nan(length(all_strains), 2);

[~, a, b] = intersect(all_strains, hit_strains);
all_data(a, 1) = hit_data(b);

[~, a, b] = intersect(all_strains, hit_strains2);
all_data(a, 2) = hit_data2(b);

%% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5183; 5184];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
gonzalez_valero_2016.orfs = all_strains;
gonzalez_valero_2016.ph = hit_data_names;
gonzalez_valero_2016.data = all_data;
gonzalez_valero_2016.dataset_ids = hit_data_ids;

%% Save

save('./gonzalez_valero_2016.mat','gonzalez_valero_2016');

%% Print out

fid = fopen('./gonzalez_valero_2016.txt','w');
write_matrix_file(fid, gonzalez_valero_2016.orfs, gonzalez_valero_2016.ph, gonzalez_valero_2016.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(gonzalez_valero_2016)
end

end

