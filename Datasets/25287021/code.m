%% Pereira~Domingues, 2014
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
pereira_domingues_2014.pmid = 25287021;

hit_data_ids = [751; 752];

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(pereira_domingues_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% All Strains
[FILENAMES{end+1}, all_orfs] = read_data('xlsread','./raw_data/chemogenomics.xlsx');

% Get the list of orfs
tested_orfs = all_orfs(2:end,1);

% Clean them up
tested_orfs = clean_orf(tested_orfs);

% Find anything that doesn't look like an ORF and remove it
inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(tested_orfs);

%% Hit Strains
[FILENAMES{end+1}, wsh] = read_data('xlsread','./raw_data/10295_2014_1519_MOESM1_ESM.xlsx');
[FILENAMES{end+1}, sh] = read_data('xlsread','./raw_data/10295_2014_1519_MOESM2_ESM.xlsx');

% Get the list of ORFs
ind = find(cellfun(@ischar, sh(10:end,11)));
hit_sh = sh(10:end, 1);
sh_data = sh(10:end, 11);
hit_sh = hit_sh(ind, 1);
sh_data = sh_data(ind, 1);

ind = find(cellfun(@ischar, wsh(10:end,11)));
hit_wsh = wsh(10:end, 1);
wsh_data = wsh(10:end, 11);
hit_wsh = hit_wsh(ind, 1);
wsh_data = wsh_data(ind, 1);

% Translate genenames to orfs
hit_wsh = clean_genename(hit_wsh);
hit_sh = clean_genename(hit_sh);
hit_wsh = translate(hit_wsh);
hit_sh = translate(hit_sh);

% Find anything that doesn't look like an ORF and remove it
inds = find(~is_orf(hit_wsh));
hit_wsh(inds) = [];
wsh_data(inds) = [];
inds = find(~is_orf(hit_sh));
hit_sh(inds) = [];
sh_data(inds) = [];

%% Data
% Make an array of zeros
hit_data_wsh = zeros(length(hit_wsh), 1)-1;
hit_data_wsh(find(~cellfun(@isempty, strfind(wsh_data, '--')))) = -2;
hit_data_sh = zeros(length(hit_sh), 1)-1;
hit_data_sh(find(~cellfun(@isempty, strfind(sh_data, '--')))) = -2;

% Match the hit names with all the names
final_data = zeros(length(tested_orfs), 2);
[~,ind1,ind2] = intersect(tested_orfs, hit_wsh);
final_data(ind1,1) = hit_data_wsh(ind2);
[~,ind1,ind2] = intersect(tested_orfs, hit_sh);
final_data(ind1,2) = hit_data_sh(ind2);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

pereira_domingues_2014.orfs = tested_orfs;
pereira_domingues_2014.ph = hit_data_names;
pereira_domingues_2014.data = final_data;
pereira_domingues_2014.dataset_ids = hit_data_ids;

%% Save

save('./pereira_domingues_2014.mat','pereira_domingues_2014');

fid = fopen('./pereira_domingues_2014.txt','w');
write_matrix_file(fid, pereira_domingues_2014.orfs, pereira_domingues_2014.ph, pereira_domingues_2014.data);
fclose(fid);

end


