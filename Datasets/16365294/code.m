%% Ohya~Morishita, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ohya_morishita_2005.pmid = 16365294;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ohya_morishita_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

[FILENAMES{end+1}, data] = read_data('importdata','./raw_data/mutant_analysis_2011_10_20.tab');
data.orfs = data.textdata(2:end,1);
data.params = data.textdata(1,2:end)';

% Get the IDs of the parameters to be retained
phenotype_name = cell(length(datasets.id),1);
phenotype_id = cell(length(datasets.id),1);
for i = 1 : length(datasets.id)
    tmp = regexp(datasets.standard_name{i}, '\|', 'split');
    t = strtrim(tmp{2});
    tmp2 = regexp(t, '\(([A-Z0-9_\-]*?)\)$', 'match');
    tmp2 = tmp2{1}(2:end-1);
    phenotype_name{i} = t;
    phenotype_id{i} = tmp2;
end

% First, get just the parameters that we want to include (exclude the CV
% parameters and the ones that are not easily interpretable)
[~,ind1,ind2] = intersect(phenotype_id, data.params);
data.params = data.params(ind2);
data.data = data.data(:,ind2);

phenotype_id = phenotype_id(ind1);
phenotype_name = phenotype_name(ind1);
hit_data_ids = datasets.id(ind1);

hit_orfs = data.orfs;
hit_data = data.data;

% Eliminate all white spaces & capitalize
hit_orfs = clean_orf(hit_orfs);

hit_orfs = translate(hit_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_orfs));
disp(hit_orfs(inds));  

% If the same strain is present more than once, average its values
[hit_orfs, hit_data] = grpstats(hit_data, hit_orfs, {'gname','mean'});

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
ohya_morishita_2005.orfs = hit_orfs;
ohya_morishita_2005.ph = hit_data_names;
ohya_morishita_2005.data = hit_data;
ohya_morishita_2005.dataset_ids = hit_data_ids;

%% Save

save('./ohya_morishita_2005.mat','ohya_morishita_2005');

%% Print out

fid = fopen('./ohya_morishita_2005.txt','w');
write_matrix_file(fid, ohya_morishita_2005.orfs, ohya_morishita_2005.ph, ohya_morishita_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ohya_morishita_2005)
end

end

