%% Fei~Yang, 2008
function FILENAMES = code()

FILENAMES = {};
fei_yang_2008.pmid = 18250201;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(fei_yang_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};


%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('textread','./raw_data/hits.txt','%s %d');

% Get the list of ORFs and the correponding data 
hit_strains = data{1};

% Get the data itself
hit_data = data{2}; % if the dataset is discrete or binary

% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [168];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
fei_yang_2008.orfs = hit_strains;
fei_yang_2008.ph = hit_data_names;
fei_yang_2008.data = hit_data;
fei_yang_2008.dataset_ids = hit_data_ids;

%% Save

save('./fei_yang_2008.mat','fei_yang_2008');

%% Print out

fid = fopen('./fei_yang_2008.txt','w');
write_matrix_file(fid, fei_yang_2008.orfs, fei_yang_2008.ph, fei_yang_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(fei_yang_2008)
end

end

