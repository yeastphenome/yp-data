%% Nislow~Hammond, 2015
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
nislow_hammond_2015.pmid = 25667933;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(nislow_hammond_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

files = {'./raw_data/TableS5_vsT1-hom1-dropsToBg.xlsx'; ...
    './raw_data/TableS6_vsT1-hom1-NaCl-dropsToBg.xlsx'; ...
    './raw_data/TableS8_vsT1-het1-dropsToBg.xlsx'; ...
    './raw_data/TableS9_vsT1-het1-NaCl-dropsToBg.xlsx'};

sheets = {'ground'; 'flight'};

for f = 1 : length(files)
    for s = 1 : length(sheets)
        
        i = (f-1)*length(sheets)+s;
        [FILENAMES{end+1}, data] = read_data('xlsread',files{f}, sheets{s});

        % Get the list of ORFs and the correponding data 
        % (this part usually changes significantly based on the format of the raw data file)
        hit_strains{i} = data(2:end,2);
        
        column_names = data(1,:);
        inds_columns = find(strncmp('Log2ratio', column_names, length('Log2ratio')));
        
        % Get the data itself
        hit_data{i} = data(2:end,inds_columns);
   
        % Eliminate all white spaces & capitalize
        hit_strains{i} = clean_orf(hit_strains{i});

        % Find anything that doesn't look like an ORF
        inds = find(~is_orf(hit_strains{i}));
        hit_strains{i}(inds) = [];
        hit_data{i}(inds, :) = [];

        hit_data{i} = cell2mat(hit_data{i});

        % If the same strain is present more than once, average its values
        [hit_strains{i}, hit_data{i}] = grpstats(hit_data{i}, hit_strains{i}, {'gname','mean'});
        
    end
end


%% Combine the data

hit_strains_all = unique(vertcat(hit_strains{:}));

n = sum(cellfun(@(x) size(x,2), hit_data));

hit_data_all = nan(length(hit_strains_all), n);
c = 1;
for i = 1 : length(hit_data)
    [~,ind1,ind2] = intersect(hit_strains{i}, hit_strains_all);
    hit_data_all(ind2,c:c+size(hit_data{i},2)-1) = hit_data{i}(ind1,:);
    c = c+size(hit_data{i},2);
end

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5276; 5277; 5278; 5279; 5280; 5281; 5282; 5283; 5284; 5285; 5286; 5287; 5288; 5289];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
nislow_hammond_2015.orfs = hit_strains_all;
nislow_hammond_2015.ph = hit_data_names;
nislow_hammond_2015.data = -hit_data_all;   % Data represent log(gen1/gen2) such that, originally, higher values correspond to more extreme depletion
nislow_hammond_2015.dataset_ids = hit_data_ids;

%% Save

save('./nislow_hammond_2015.mat','nislow_hammond_2015');

%% Print out

fid = fopen('./nislow_hammond_2015.txt','w');
write_matrix_file(fid, nislow_hammond_2015.orfs, nislow_hammond_2015.ph, nislow_hammond_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(nislow_hammond_2015)
end

end

