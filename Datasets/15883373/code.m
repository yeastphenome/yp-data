%% Xie~Huang, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
xie_huang_2005.pmid = 15883373;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(xie_huang_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};


%% 
[FILENAMES{end+1}, data1.raw] = read_data('xlsread','./raw_data/rapa_cellarray_all.xls', 'Rapa10nm');
[FILENAMES{end+1}, data2.raw] = read_data('xlsread','./raw_data/rapa_cellarray_all.xls', 'Rapa30nm');

% Get indices of the data columns
ind_data1 = find(strcmp('Rapa/DMSO', data1.raw(1,:)));
ind_data2 = find(strcmp('Rapa/DMSO', data2.raw(1,:)));

hit_orfs1 = data1.raw(:,3);
hit_orfs2 = data2.raw(:,3);

hit_data1 = data1.raw(:,ind_data1);
hit_data2 = data2.raw(:,ind_data2);

% Eliminate anything that doesn't look like an ORF
hit_orfs1 = clean_orf(hit_orfs1);
hit_orfs2 = clean_orf(hit_orfs2);

inds = find(~is_orf(hit_orfs1));
disp(hit_orfs1(inds));
hit_orfs1(inds) = [];
hit_data1(inds,:) = [];

inds = find(~is_orf(hit_orfs2));
disp(hit_orfs2(inds));
hit_orfs2(inds) = [];
hit_data2(inds,:) = [];

hit_data1 = cell2mat(hit_data1);
hit_data2 = cell2mat(hit_data2);

% Average data for identical ORFs that appear multiple times
[hit_orfs1,hit_data1] = grpstats(hit_data1, hit_orfs1, {'gname','mean'});
[hit_orfs2,hit_data2] = grpstats(hit_data2, hit_orfs2, {'gname','mean'});

hit_orfs = unique([hit_orfs1; hit_orfs2]);
hit_data = nan(length(hit_orfs),2);

[~,ind1,ind2] = intersect(hit_orfs1, hit_orfs);
hit_data(ind2,1) = hit_data1(ind1);
[~,ind1,ind2] = intersect(hit_orfs2, hit_orfs);
hit_data(ind2,2) = hit_data2(ind1);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [55; 56];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

xie_huang_2005.orfs = hit_orfs;
xie_huang_2005.data = hit_data;
xie_huang_2005.ph = hit_data_names;
xie_huang_2005.dataset_ids = hit_data_ids;

%% Save

save('./xie_huang_2005.mat','xie_huang_2005');


%% Print out

fid = fopen('./xie_huang_2005.txt','w');
write_matrix_file(fid, xie_huang_2005.orfs, xie_huang_2005.ph, xie_huang_2005.data);
fclose(fid);

end
