%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
fujita_iwahashi_2006.pmid = 16879425;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(fujita_iwahashi_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

for i = 1 : 3
    [FILENAMES{end+1}, data{i}] = read_data('readtable',['./raw_data/Table' num2str(i) '.txt'], 'HeaderLines', 1);

    % Get the list of ORFs and the correponding data 
    % (this part usually changes significantly based on the format of the raw data file)
    hit_strains{i} = data{i}.ORF;

    % Eliminate all white spaces & capitalize
    hit_strains{i} = clean_orf(hit_strains{i});

    % If in gene name form, transform into ORF name
    hit_strains{i} = translate(hit_strains{i});

    % Find anything that doesn't look like an ORF
    inds = find(~is_orf(hit_strains{i}));
    disp(hit_strains{i}(inds));  
    hit_strains{i}(inds) = [];

end

hit_strains_all = unique(vertcat(hit_strains{:}));
hit_data_all = zeros(length(hit_strains_all),3);

for i = 1 : 3
    [~,ind1,ind2] = intersect(hit_strains_all, hit_strains{i});
    hit_data_all(ind1,i) = -1;
end

% If the same strain is present more than once, average its values
[hit_strains_all, hit_data_all] = grpstats(hit_data_all, hit_strains_all, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [497 5004 5005]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
fujita_iwahashi_2006.orfs = hit_strains_all;
fujita_iwahashi_2006.ph = hit_data_names;
fujita_iwahashi_2006.data = hit_data_all;
fujita_iwahashi_2006.dataset_ids = hit_data_ids;

%% Save

save('./fujita_iwahashi_2006.mat','fujita_iwahashi_2006');

%% Print out

fid = fopen('./fujita_iwahashi_2006.txt','w');
write_matrix_file(fid, fujita_iwahashi_2006.orfs, fujita_iwahashi_2006.ph, fujita_iwahashi_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(fujita_iwahashi_2006)
end

end

