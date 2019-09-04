%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
pautasso_rossi_2016.pmid = 27188886;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(pautasso_rossi_2016.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Supplementary_Table_1.xlsx', 'ORFs 2&lt;Z&lt;2');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)

hit_strains = {};
hit_data = {};

for g = 1:4
    
    i = g*2-1;
    
    hit_strains_this = [data(5:end,i); data(5:end,i+1)];
    hit_data_this = [-ones(size(data(5:end,i))); ones(size(data(5:end,i+1)))];

    inds = find(~cellfun(@isstr, hit_strains_this));
    hit_strains_this(inds) = [];
    hit_data_this(inds) = [];

    hit_strains_this = clean_genename(hit_strains_this);
    hit_strains_this = translate(hit_strains_this);
    inds = find(~is_orf(hit_strains_this));
    disp(hit_strains_this(inds)); 
    
    hit_strains_this(find(strcmp('FLO8', hit_strains_this))) = {'YER109C'};
    
    [hit_strains_this, hit_data_this] = grpstats(hit_data_this, hit_strains_this, {'gname','mean'});

    hit_strains{g} = hit_strains_this;
    hit_data{g} = hit_data_this;
end
   
hit_strains_all = unique(vertcat(hit_strains{:}));
hit_data_all = zeros(length(hit_strains_all),4);

for g = 1:4
    
    [~,ind1,ind2] = intersect(hit_strains_all, hit_strains{g});
    hit_data_all(ind1,g) = hit_data{g}(ind2);
end

% If the same strain is present more than once, average its values

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16247; 16248; 16249; 16250];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
pautasso_rossi_2016.orfs = hit_strains_all;
pautasso_rossi_2016.ph = hit_data_names;
pautasso_rossi_2016.data = hit_data_all;
pautasso_rossi_2016.dataset_ids = hit_data_ids;

%% Save

save('./pautasso_rossi_2016.mat','pautasso_rossi_2016');

%% Print out

fid = fopen('./pautasso_rossi_2016.txt','w');
write_matrix_file(fid, pautasso_rossi_2016.orfs, pautasso_rossi_2016.ph, pautasso_rossi_2016.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(pautasso_rossi_2016)
end

end

