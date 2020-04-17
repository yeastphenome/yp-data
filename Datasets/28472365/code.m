%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
maclean_zhang_2017.pmid = 28472365;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(maclean_zhang_2017.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

hit_data_ids = [16395, 16397, 16399, 16396, 16400, 16398]';
sheets = {'High Temp. (40C)','EtOH','H2O2','NaCl','CoCL2','SO'};

hit_strains_all = {};
hit_data_all = {};

for i = 1:length(sheets)

    [FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supplementary_data_6.xlsx', sheets{i});

    % Get the list of ORFs and the correponding data 
    % (this part usually changes significantly based on the format of the raw data file)
    hit_strains = data(2:end,1);

    % Remove the underscore annotations
    t = regexp(hit_strains, '_', 'split');
    hit_strains = cellfun(@(v)v(1), t);

    % Get the data itself
    hit_data = data(2:end,3);

    % Eliminate all white spaces & capitalize
    hit_strains = clean_genename(hit_strains);

    % If in gene name form, transform into ORF name
    hit_strains = translate(hit_strains);

    % Find anything that doesn't look like an ORF
    inds = find(~is_orf(hit_strains));
    disp(hit_strains(inds));  

    hit_data = cell2mat(hit_data);

    % If the same strain is present more than once, average its values
    [hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});
    
    hit_strains_all{i} = hit_strains;
    hit_data_all{i} = hit_data;
    
end

hit_strains = unique(vertcat(hit_strains_all{:}));
hit_data = zeros(length(hit_strains), length(sheets));

for i = 1 : length(sheets)
    [~,ind1,ind2] = intersect(hit_strains, hit_strains_all{i});
    hit_data(ind1,i) = hit_data_all{i}(ind2,:);
end


%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
maclean_zhang_2017.orfs = hit_strains;
maclean_zhang_2017.ph = hit_data_names;
maclean_zhang_2017.data = hit_data;
maclean_zhang_2017.dataset_ids = hit_data_ids;

%% Save

save('./maclean_zhang_2017.mat','maclean_zhang_2017');

%% Print out

fid = fopen('./maclean_zhang_2017.txt','w');
write_matrix_file(fid, maclean_zhang_2017.orfs, maclean_zhang_2017.ph, maclean_zhang_2017.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(maclean_zhang_2017)
end

end

