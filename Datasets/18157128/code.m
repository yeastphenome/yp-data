%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
delneri_oliver_2008.pmid = 18157128;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(delneri_oliver_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table 4s_NewForPublishing.xlsx');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
ind_strains = find(strcmp('ORF', data(11,:)));
for i = 1:length(ind_strains)/2
    c = i*2-1;
    hit_strains{i} = [data(13:end,ind_strains(c)); data(13:end, ind_strains(c+1))];
end

% Get the data itself
ind_data = ind_strains+1;
for i = 1:length(ind_data)/2
    c = i*2-1;
    hit_data{i} = [data(13:end,ind_data(c)); data(13:end, ind_data(c+1))];
end
   
% Eliminate all white spaces & capitalize
for i = 1 : 4
    
    inds = find(cellfun(@isnumeric, hit_strains{i}));
    hit_strains{i}(inds) = [];
    hit_data{i}(inds) = [];
    
    inds = find(~cellfun(@isnumeric, hit_data{i}));
    hit_strains{i}(inds) = [];
    hit_data{i}(inds) = [];
    
    hit_data{i} = cell2mat(hit_data{i});
    
    hit_strains{i} = clean_orf(hit_strains{i});

    % If in gene name form, transform into ORF name
    hit_strains{i} = translate(hit_strains{i});

    % Find anything that doesn't look like an ORF
    inds = find(~is_orf(hit_strains{i}));
    disp(hit_strains{i}(inds));  
    
end


hit_strains_all = unique(cat(1,hit_strains{:}));
hit_data_all = nan(length(hit_strains_all),4);

for i = 1 : 4
    
    % If the same strain is present more than once, average its values
    [t1, t2] = grpstats(hit_data{i}, hit_strains{i}, {'gname','mean'});
    
    [~,ind1,ind2] = intersect(hit_strains_all, t1);
    hit_data_all(ind1,i) = t2(ind2);
end
    
% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11813, 11815, 11816, 11814]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
delneri_oliver_2008.orfs = hit_strains_all;
delneri_oliver_2008.ph = hit_data_names;
delneri_oliver_2008.data = hit_data_all;
delneri_oliver_2008.dataset_ids = hit_data_ids;

%% Save

save('./delneri_oliver_2008.mat','delneri_oliver_2008');

%% Print out

fid = fopen('./delneri_oliver_2008.txt','w');
write_matrix_file(fid, delneri_oliver_2008.orfs, delneri_oliver_2008.ph, delneri_oliver_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(delneri_oliver_2008)
end

end

