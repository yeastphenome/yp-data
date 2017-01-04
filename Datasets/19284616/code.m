%% Thorsen~Tamas, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
thorsen_tamas_2009.pmid = 19284616;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(thorsen_tamas_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data
% As (III) - 0.5 mM & 24 hr

sheets = {'0.5 As 24','1 As 24h','1.5 As 24h','0.5As 48h','1 As 48h','1.5 As 48h'};

for s = 1 : length(sheets)
    
    [FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/As(III) list avg.xlsx', sheets{s});

    % Get the list of ORFs and the correponding data 
    hit_strains{s} = data(2:end,11);

    % Get the data itself
    hit_data{s} = cell2mat(data(2:end,10));

    % Eliminate all white spaces & capitalize
    hit_strains{s} = clean_orf(hit_strains{s});

    % Find anything that doesn't look like an ORF
    inds = find(~is_orf(hit_strains{s}));
    hit_strains{s}(inds) = [];
    hit_data{s}(inds) = [];

    % If the same strain is present more than once, average its values
    [hit_strains{s}, hit_data{s}] = grpstats(hit_data{s}, hit_strains{s}, {'gname','mean'});

end


% Cd
sheets = {'Cd 75 24h','Cd 100 24h','Cd 150 24h','Cd 75 48h','Cd 100 48h','Cd 150 48h'};

for s = 1 : length(sheets)
    
    [FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Cd list avg.xlsx', sheets{s});

    % Get the list of ORFs and the correponding data 
    hit_strains2{s} = data(2:end,10);

    % Get the data itself
    hit_data2{s} = data(2:end,4);
    hit_data2{s}(~cellfun(@(x) any(isnumeric(x(:))), hit_data2{s})) = {NaN};
    hit_data2{s} = cell2mat(hit_data2{s});

    % Eliminate all white spaces & capitalize
    hit_strains2{s} = clean_orf(hit_strains2{s});

    % Find anything that doesn't look like an ORF
    inds = find(~is_orf(hit_strains2{s}));
    hit_strains2{s}(inds) = [];
    hit_data2{s}(inds) = [];

    % If the same strain is present more than once, average its values
    [hit_strains2{s}, hit_data2{s}] = grpstats(hit_data2{s}, hit_strains2{s}, {'gname','mean'});
    
end

%% Combine all datasets
hit_strains_all = unique([vertcat(hit_strains{:}); vertcat(hit_strains2{:})]);
hit_data_all = nan(length(hit_strains_all), 12);

for s = 1 : 6
    [~, ind1, ind2] = intersect(hit_strains_all, hit_strains{s});
    hit_data_all(ind1, s) = hit_data{s}(ind2);
    
    [~, ind1, ind2] = intersect(hit_strains_all, hit_strains2{s});
    hit_data_all(ind1, s+6) = hit_data2{s}(ind2);
end

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [1320; 5362; 5363; 5364; 5365; 5366; 1319; 5367; 5368; 5369; 5370; 5371];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
thorsen_tamas_2009.orfs = hit_strains_all;
thorsen_tamas_2009.ph = hit_data_names;
thorsen_tamas_2009.data = hit_data_all;
thorsen_tamas_2009.dataset_ids = hit_data_ids;

%% Save

save('./thorsen_tamas_2009.mat','thorsen_tamas_2009');

%% Print out

fid = fopen('./thorsen_tamas_2009.txt','w');
write_matrix_file(fid, thorsen_tamas_2009.orfs, thorsen_tamas_2009.ph, thorsen_tamas_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(thorsen_tamas_2009)
end

end