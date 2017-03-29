%% Gaupel~Tenniswood, 2014
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
gaupel_tenniswood_2004.pmid = 24968945;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(gaupel_tenniswood_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data_sens.raw] = read_data('xlsread','./raw_data/1471-2164-15-528-s1.xlsx', 'Top sensitive strains');

hits_orfs1 = data_sens.raw(2:end,1);
hits_data1 = -cell2mat(data_sens.raw(2:end,3));

hits_orfs1 = clean_orf(hits_orfs1);

hits_orfs1(ismember(hits_orfs1, {'YJL206-A'})) = {'YJL206C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs1));
disp(hits_orfs1(inds));  

%% Load data (2)

[FILENAMES{end+1}, data_res.raw] = read_data('xlsread','./raw_data/1471-2164-15-528-s1.xlsx', 'Top resistant strains');

hits_orfs2 = data_res.raw(2:end,1);
hits_data2 = -cell2mat(data_res.raw(2:end,3));

hits_orfs2 = clean_orf(hits_orfs2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs2));
disp(hits_orfs2(inds));  

%% Combine part1 and part2

hits_orfs = [hits_orfs1; hits_orfs2];
hits_data = [hits_data1; hits_data2];

[hits_orfs, hits_data] = grpstats(hits_data, hits_orfs,{'gname','mean'});

%% Tested strains

% Load tested strains
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/CompleteDeletionLibrary.xlsx');
tested_orfs = tested.raw(2:end,1);

% Eliminate all white spaces & capitalize
tested_orfs = clean_orf(tested_orfs);

% If possible, fix the typo
tested_orfs(ismember(tested_orfs, {'YAR002AW'})) = {'YAR002W-A'};
tested_orfs(ismember(tested_orfs, {'YOLO57W'})) = {'YOL057W'};
tested_orfs(ismember(tested_orfs, {'YKLO72W'})) = {'YKL072W'};
tested_orfs(ismember(tested_orfs, {'YJL206-A'})) = {'YJL206C-A'};
tested_orfs(ismember(tested_orfs, {'YLR287-A'})) = {'YLR287C-A'};
tested_orfs(ismember(tested_orfs, {'YFL033AC'})) = {'YFL033C-A'};
tested_orfs(ismember(tested_orfs, {'YOLO62C'})) = {'YOL062C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds)); 

% Take the unique set
tested_orfs = unique(tested_orfs);

% Make sure that all the hits are part of the tested set
[missing, ix] = setdiff(hits_orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [569];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
gaupel_tenniswood_2004.orfs = tested_orfs;
gaupel_tenniswood_2004.ph = hit_data_names;
gaupel_tenniswood_2004.data = zeros(length(gaupel_tenniswood_2004.orfs),length(gaupel_tenniswood_2004.ph));
gaupel_tenniswood_2004.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
gaupel_tenniswood_2004.data(ind2,:) = hits_data(ind1,:);

%% Save

save('./gaupel_tenniswood_2004.mat','gaupel_tenniswood_2004');

%% Print out

fid = fopen('./gaupel_tenniswood_2004.txt','w');
write_matrix_file(fid, gaupel_tenniswood_2004.orfs, gaupel_tenniswood_2004.ph, gaupel_tenniswood_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(gaupel_tenniswood_2004)
end

end
