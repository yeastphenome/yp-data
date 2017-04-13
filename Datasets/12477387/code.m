%% Zhang~Schneider, 2002
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
zhang_schneider_2002.pmid = 12477387;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(zhang_schneider_2002.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};


%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Zhang et al supplemental data.xlsx', 'Size Data');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);

% Get the data itself
hit_data = data(:,5);   % keep just the mean
    
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% A few corrections
hit_strains(strcmp('TAL004W', hit_strains)) = {'YAL004W'};
hit_strains(strcmp('YELOO1C', hit_strains)) = {'YEL001C'};
hit_strains(strcmp('KL187C', hit_strains)) = {'YKL187C'};
hit_strains(strcmp('YMR41W', hit_strains)) = {'YMR241W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If not possible, eliminate the entry
hit_strains(inds) = [];
hit_data(inds,:) = [];

hit_data = cell2mat(hit_data);

hit_data2 = nan(length(hit_strains),2);

% To separate HOM from HET data, remove all genes that don't appear on
% today's HOM collection from Open Biosystems
[FILENAMES{end+1}, hom] = read_data('xlsread','./extras/Homozygous_diploid_obs_v7.0.xlsx', 'DATA');
hom_orfs = unique(translate(clean_orf(hom(2:end,2))));

inds = find(ismember(hit_strains, hom_orfs));
hit_data2(inds,1) = hit_data(inds,1);
inds = find(~ismember(hit_strains, hom_orfs));
hit_data2(inds,2) = hit_data(inds,1);

hit_data = hit_data2;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [477 5384]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
zhang_schneider_2002.orfs = hit_strains;
zhang_schneider_2002.ph = hit_data_names;
zhang_schneider_2002.data = hit_data;
zhang_schneider_2002.dataset_ids = hit_data_ids;

%% Save

save('./zhang_schneider_2002.mat','zhang_schneider_2002');

%% Print out

fid = fopen('./zhang_schneider_2002.txt','w');
write_matrix_file(fid, zhang_schneider_2002.orfs, zhang_schneider_2002.ph, zhang_schneider_2002.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(zhang_schneider_2002)
end

end

