%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hoepfner_movva_2014.pmid = 24360837;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hoepfner_movva_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('read_matrix_file','./raw_data/HOP_scores.txt');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data.labels_row;

% Get the data itself
hit_data = data.data; % if the dataset is discrete or binary
  
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

hit_strains(ismember(hit_strains, {'YBR160WAS'})) = {'YBR160W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

% HET data
[FILENAMES{end+1}, data_hip] = read_data('read_matrix_file','./raw_data/HIP_scores.txt');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains_hip = data_hip.labels_row;

% Get the data itself
hit_data_hip = data_hip.data; % if the dataset is discrete or binary
  
% Eliminate all white spaces & capitalize
hit_strains_hip = clean_orf(hit_strains_hip);

hit_strains_hip(ismember(hit_strains_hip, {'YBR160WAS'})) = {'YBR160W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_hip));
disp(hit_strains_hip(inds));  

hit_strains_hip(inds) = [];
hit_data_hip(inds,:) = [];

%% Map dataset IDs to data columns

% First, eliminate all z-scores
inds = find(~cellfun(@isempty, regexp(data.labels_col, 'z-score')));
data.labels_col(inds) = [];
hit_data(:,inds) = [];

inds = find(~cellfun(@isempty, regexp(data_hip.labels_col, 'z-score')));
data_hip.labels_col(inds) = [];
hit_data_hip(:,inds) = [];

str = '"%s scores for Exp. %d_%f_HOP_%s"';

type = cell(size(data.labels_col));
exp = nan(size(data.labels_col));
conc = nan(size(data.labels_col));
etc = cell(size(data.labels_col));
for i = 1 : length(data.labels_col)
    C = textscan(data.labels_col{i}, str);
    type(i) = C{1};
    exp(i) = C{2};
    conc(i) = C{3};
    etc(i) = C{4};
end

type_hip = cell(size(data_hip.labels_col));
exp_hip = nan(size(data_hip.labels_col));
conc_hip = nan(size(data_hip.labels_col));
etc_hip = cell(size(data_hip.labels_col));
for i = 1 : length(data_hip.labels_col)
    C = textscan(data_hip.labels_col{i}, str);
    type_hip(i) = C{1};
    exp_hip(i) = C{2};
    conc_hip(i) = C{3};
    if ~isempty(C{4})
        etc_hip(i) = C{4};
    end
end

% Load CMD to dataset mapping for the known compounds
T = readtable('./extras/type_cmb_dose_dataset.txt','Delimiter','\t');
T_hip = T;

% Only focus on the known compounds
inds = find(~ismember(exp, T.CMB));
type(inds) = [];
exp(inds) = [];
conc(inds) = [];
etc(inds) = [];
hit_data(:,inds) = [];

hit_data_ids = nan(size(type));
for i = 1 : length(exp)
    inds = find(strcmp(type{i}, T.Type) & T.CMB==exp(i) & T.Dose==conc(i));
    hit_data_ids(i) = unique(T.DatasetHOP(inds));
end


% Only focus on the known compounds
inds = find(~ismember(exp_hip, T_hip.CMB));
type_hip(inds) = [];
exp_hip(inds) = [];
conc_hip(inds) = [];
etc_hip(inds) = [];
hit_data_hip(:,inds) = [];

hit_data_ids_hip = nan(size(type_hip));
for i = 1 : length(exp_hip)
    inds = find(strcmp(type_hip{i}, T_hip.Type) & T_hip.CMB==exp_hip(i) & T_hip.Dose==conc_hip(i));
    if ~isempty(inds)
        hit_data_ids_hip(i) = unique(T.DatasetHIP(inds));
    end
end

inds = find(isnan(hit_data_ids_hip));
hit_data_ids_hip(inds) = [];
hit_data_hip(:,inds) = [];


%%

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% If the same dataset is present more than once, average its values (as it
% represents replicates)
[hit_data_ids, hit_data] = grpstats(hit_data', hit_data_ids, {'gname','mean'});
hit_data = hit_data';
hit_data_ids = cellfun(@str2num, hit_data_ids);


% If the same strain is present more than once, average its values
[hit_strains_hip, hit_data_hip] = grpstats(hit_data_hip, hit_strains_hip, {'gname','mean'});

% If the same dataset is present more than once, average its values (as it
% represents replicates)
[hit_data_ids_hip, hit_data_hip] = grpstats(hit_data_hip', hit_data_ids_hip, {'gname','mean'});
hit_data_hip = hit_data_hip';
hit_data_ids_hip = cellfun(@str2num, hit_data_ids_hip);

% Combine the HIP and HOP data together

hit_strains_all = unique([hit_strains; hit_strains_hip]);
hit_data_all = nan(length(hit_strains_all), length(hit_data_ids)+length(hit_data_ids_hip));

[~,ind1,ind2] = intersect(hit_strains_all, hit_strains);
hit_data_all(ind1,1:length(hit_data_ids)) = hit_data(ind2,:);

[~,ind1,ind2] = intersect(hit_strains_all, hit_strains_hip);
hit_data_all(ind1,length(hit_data_ids)+1:end) = hit_data_hip(ind2,:);

hit_data_ids_all = [hit_data_ids; hit_data_ids_hip];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids_all);
hit_data_names = cell(size(hit_data_ids_all));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
hoepfner_movva_2014.orfs = hit_strains_all;
hoepfner_movva_2014.ph = hit_data_names;
hoepfner_movva_2014.data = hit_data_all;
hoepfner_movva_2014.dataset_ids = hit_data_ids_all;

% Eliminate all ORFs that have very few values (are likely to be essential -- in HOM data it's hard to exclude them otherwise)
% numvalues = sum(~isnan(hoepfner_movva_2014.data),2);
% inds = find(numvalues<round(0.25 * length(hoepfner_movva_2014.ph)));
% hoepfner_movva_2014.orfs(inds) = [];
% hoepfner_movva_2014.data(inds,:) = [];

%% Save

save('./hoepfner_movva_2014.mat','hoepfner_movva_2014');

%% Print out

fid = fopen('./hoepfner_movva_2014.txt','w');
write_matrix_file(fid, hoepfner_movva_2014.orfs, hoepfner_movva_2014.ph, hoepfner_movva_2014.data);
fclose(fid);

%% Save to DB (admin)

% addpath(genpath('../../Private-Utils/'));
% if exist('save_data_to_db.m')
%     res = save_data_to_db(hoepfner_movva_2014)
% end

% Special output for DB (given the size)
fid = fopen('db.txt','w');
for i = 1 : length(hoepfner_movva_2014.dataset_ids)
    for j = 1 : length(hoepfner_movva_2014.orfs)
        if ~isnan(hoepfner_movva_2014.data(j,i))
            fprintf(fid, '%s\t%d\t%.3f\n', hoepfner_movva_2014.orfs{j}, hoepfner_movva_2014.dataset_ids(i), hoepfner_movva_2014.data(j,i));
        end
    end
    i
end
fclose(fid);


end

