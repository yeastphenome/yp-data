%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
deutschbauer_davis_2002.pmid = 12432101;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(deutschbauer_davis_2002.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

% Load quant for hits data
[FILENAMES{end+1}, data_def] = read_data('readtable','./raw_data/sporulation_deficient.txt');
[FILENAMES{end+1}, data_prof] = read_data('readtable','./raw_data/sporulation_proficient.txt');
[FILENAMES{end+1}, data_germ] = read_data('readtable','./raw_data/germination.txt');

sporulation_strains = [data_def.Orf; data_prof.Orf];
sporulation_data = [data_def{:,{'Spo1','Spo2'}}; data_prof{:,{'Spo1','Spo2'}}];

% Taking the reciprocal so that higher values = higher growth
sporulation_data = 1./sporulation_data;

% Average the 2 experiments
sporulation_data = nanmean(sporulation_data, 2);

germ_strains = data_germ.Orf;
germ_data = nanmean(data_germ{:,{'Germ1','Germ2'}},2);

hit_strains = unique([sporulation_strains; germ_strains]);
hit_data = nan(length(hit_strains),2);
[~,ind1,ind2] = intersect(hit_strains, sporulation_strains);
hit_data(ind1,1) = sporulation_data(ind2);
[~,ind1,ind2] = intersect(hit_strains, germ_strains);
hit_data(ind1,2) = germ_data(ind2);

  
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

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
hit_data_ids = [478; 16004];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('readtable','./raw_data/SpoGerm_RawData.txt');

tested_strains = tested_strains.ORF;

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If in gene name form, transform into ORF name
tested_strains = translate(tested_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Separate essentials and non-essentials
load essential_genes_151215.mat
inds_essential = find(ismember(tested_strains, essential_genes));
tested_strains(inds_essential) = [];

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
deutschbauer_davis_2002.orfs = tested_strains;
deutschbauer_davis_2002.ph = hit_data_names;
deutschbauer_davis_2002.data = ones(length(deutschbauer_davis_2002.orfs),length(deutschbauer_davis_2002.ph));   % NOTE: WT = 1
deutschbauer_davis_2002.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, deutschbauer_davis_2002.orfs);
deutschbauer_davis_2002.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./deutschbauer_davis_2002.mat','deutschbauer_davis_2002');

%% Print out

fid = fopen('./deutschbauer_davis_2002.txt','w');
write_matrix_file(fid, deutschbauer_davis_2002.orfs, deutschbauer_davis_2002.ph, deutschbauer_davis_2002.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(deutschbauer_davis_2002)
end

end

