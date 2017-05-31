%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
snoek_steensma_2006.pmid = 16630279;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(snoek_steensma_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('textscan','./raw_data/hom_hits.txt', '%s %s %s', 'delimiter','\t');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data{1};
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains = unique(hit_strains);

hit_data = -ones(size(hit_strains));

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [498];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('textscan','./raw_data/Homo_diploids_041902.txt', '%s %s %*[^\n]', 'delimiter','\t');

tested_strains = tested_strains{2};

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'YELOO1C'})) = {'YEL001C'};

% If in gene name form, transform into ORF name
tested_strains = translate(tested_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

tested_strains(inds) = [];

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
snoek_steensma_2006.orfs = tested_strains;
snoek_steensma_2006.ph = hit_data_names;
snoek_steensma_2006.data = zeros(length(snoek_steensma_2006.orfs),length(snoek_steensma_2006.ph));
snoek_steensma_2006.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, snoek_steensma_2006.orfs);
snoek_steensma_2006.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./snoek_steensma_2006.mat','snoek_steensma_2006');

%% Print out

fid = fopen('./snoek_steensma_2006.txt','w');
write_matrix_file(fid, snoek_steensma_2006.orfs, snoek_steensma_2006.ph, snoek_steensma_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(snoek_steensma_2006)
end

end

