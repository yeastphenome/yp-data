%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
galardini_beltrao_2019.pmid = 31885205;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(galardini_beltrao_2019.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('readtable','./raw_data/Table_EV2_ko_scores.txt', 'Delimiter', '\t');

inds = find(strcmp('S288C', data.strain));
data = data(inds,:);

t = grpstats(data, {'condition','gene'}, {'gname','mean'}, 'DataVars', 'score');

% Get the list of ORFs and the correponding data 
hit_strains = unique(t.gene);
hit_conditions = unique(t.condition);

% Get the data itself
hit_data = nan(length(hit_strains), length(hit_conditions));

for i = 1 : length(hit_conditions)
    inds = find(strcmp(hit_conditions{i}, t.condition));
    [~,ind1,ind2] = intersect(t.gene(inds), hit_strains);
    hit_data(ind2,i) = t.mean_score(inds(ind1));
end
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YLR287-A'})) = {'YLR287C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% % Normalize by WT
% hit_data = hit_data ./ repmat(hit_data(inds,:), length(hit_strains),1);
% 
% hit_data(inds,:) = [];
% hit_strains(inds) = [];

% Load the dataset IDs
[FILENAMES{end+1}, extras] = read_data('xlsread','./raw_data/msb198831-sup-0002-tableev1.xlsx', 'Sheet1');
[~,ind1,ind2] = intersect(extras(:,1), hit_conditions);
hit_data_ids = zeros(length(hit_conditions),1);
hit_data_ids(ind2) = cell2mat(extras(ind1,2));

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
galardini_beltrao_2019.orfs = hit_strains;
galardini_beltrao_2019.ph = hit_data_names;
galardini_beltrao_2019.data = hit_data;
galardini_beltrao_2019.dataset_ids = hit_data_ids;

%% Save

save('./galardini_beltrao_2019.mat','galardini_beltrao_2019');

%% Print out

fid = fopen('./galardini_beltrao_2019.txt','w');
write_matrix_file(fid, galardini_beltrao_2019.orfs, galardini_beltrao_2019.ph, galardini_beltrao_2019.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(galardini_beltrao_2019)
end

end

