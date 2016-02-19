%% Ju~Xie, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ju_xie_2008.pmid = 18070918;

phenotypes = {'abundance of Rpn4'};
treatments = {'standard'};

% Load tested
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/mat_a_041902.xlsx', 'mat_a_041902');
tested_orfs = data.raw(3:end,2);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

% Validate tested
tested_orfs(strcmp('YLR287-A', tested_orfs)) = {'YLR287C-A'};
tested_orfs = unique(upper(strtrim(tested_orfs)));
tested_orfs(~is_orf(tested_orfs)) = [];

% Load data
fid = fopen('./raw_data/hits.txt','r');
C = textscan(fid, '%s');
fclose(fid);

genenames = C{1};

orfs = translate(genenames);
orfs(~is_orf(orfs)) = [];

data = ones(size(orfs));

missing_orfs = setdiff(orfs, tested_orfs);

ju_xie_2008.orfs = tested_orfs;
ju_xie_2008.data = zeros(length(tested_orfs),1);
[~,ind1,ind2] = intersect(tested_orfs, orfs);
ju_xie_2008.data(ind1) = data(ind2);
ju_xie_2008.ph = strcat(phenotypes, '; ', treatments);

save('./ju_xie_2008.mat','ju_xie_2008');
return;

% Save data into database
dt = ju_xie_2008;

datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% % Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

fid = fopen('./ju_xie_2008.txt','w');
write_matrix_file(fid, ju_xie_2008.orfs, ju_xie_2008.ph, ju_xie_2008.data);
fclose(fid);

end

