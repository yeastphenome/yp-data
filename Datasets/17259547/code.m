%% Zakrzewska~Klis, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
zakrzewska_klis_2007.pmid = 17259547;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(zakrzewska_klis_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

sheets = {'5h_non-essential','5h_essential','9h_non-essential','9h_essential'};

for s = 1 : length(sheets)

    [FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/SD1.xlsx', sheets{s});
    
    % Get the list of ORFs and the correponding data
    hit_strains{s} = data(2:end,1);
    
    % Get the data itself
    hit_data{s} = cell2mat(data(2:end,2));
  
    % Eliminate all white spaces & capitalize
    hit_strains{s} = clean_orf(hit_strains{s});
    
    % Find anything that doesn't look like an ORF
    inds = find(~is_orf(hit_strains{s}));
    hit_strains{s}(inds) = [];
    hit_data{s}(inds) = [];
    
    % If the same strain is present more than once, average its values
    [hit_strains{s}, hit_data{s}] = grpstats(hit_data{s}, hit_strains{s}, {'gname','mean'});

end

%% Combine all datasets
hit_strains_all = unique(vertcat(hit_strains{:}));
hit_data_all = nan(length(hit_strains_all), 4);

for s = 1 : 4
    [~, ind1, ind2] = intersect(hit_strains_all, hit_strains{s});
    hit_data_all(ind1, s) = hit_data{s}(ind2);
end

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [605; 5291; 5292; 5293];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
zakrzewska_klis_2007.orfs = hit_strains_all;
zakrzewska_klis_2007.ph = hit_data_names;
zakrzewska_klis_2007.data = hit_data_all;
zakrzewska_klis_2007.dataset_ids = hit_data_ids;

%% Save

save('./zakrzewska_klis_2007.mat','zakrzewska_klis_2007');

%% Print out

fid = fopen('./zakrzewska_klis_2007.txt','w');
write_matrix_file(fid, zakrzewska_klis_2007.orfs, zakrzewska_klis_2007.ph, zakrzewska_klis_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(zakrzewska_klis_2007)
end

end

