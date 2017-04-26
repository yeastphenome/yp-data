%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
plateau_lazard_2017.pmid = 28303947;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(plateau_lazard_2017.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/srep44761-s2.xlsx', 'data');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);

% Get the data itself
inds_data = find(ismember(data(1,:), {'Avg. 12','Avg.20'}));
hit_data = data(:,inds_data); % if the dataset is discrete or binary

% Only keep the deletions
inds = find(strcmp('DAMP', data(:,3)));
hit_strains(inds) = [];
hit_data(inds,:) = [];
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {NaN};
hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5385; 5386];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
plateau_lazard_2017.orfs = hit_strains;
plateau_lazard_2017.ph = hit_data_names;
plateau_lazard_2017.data = hit_data;
plateau_lazard_2017.dataset_ids = hit_data_ids;

%% Save

save('./plateau_lazard_2017.mat','plateau_lazard_2017');

%% Print out

fid = fopen('./plateau_lazard_2017.txt','w');
write_matrix_file(fid, plateau_lazard_2017.orfs, plateau_lazard_2017.ph, plateau_lazard_2017.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(plateau_lazard_2017)
end

end

