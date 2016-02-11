%% Suzuki~Yoshida, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
suzuki_yoshida_2011.pmid = 21601516;

phenotypes = {'glutathione abundance'};
treatments = {'standard'};

% Load data
[FILENAMES{end+1}, data.raw] = readdata('xlsread','./raw_data/hits_data.xlsx');

hits_genes = data.raw(1:end,1);
hits_data = cell2mat(data.raw(1:end,2));

hits_genes = cleanGenename(hits_genes);
hits_orfs = translate(hits_genes);

[t,t2] = grpstats(hits_data,hits_orfs,{'gname','mean'});
hits_orfs = t;
hits_data = t2;

% Load tested
[FILENAMES{end+1}, tested.raw] = readdata('xlsread','./raw_data/YKOmatalpha_GSH_list070508.xlsx');
tested_orfs = tested.raw(3:end,3);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(upper(cleanOrf(tested_orfs)));

[missing, ix] = setdiff(hits_orfs, tested_orfs);

suzuki_yoshida_2011.orfs = tested_orfs;
suzuki_yoshida_2011.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
suzuki_yoshida_2011.data(ind2,:) = hits_data(ind1,:);

suzuki_yoshida_2011.ph = strcat(phenotypes, '; ', treatments);

save('./suzuki_yoshida_2011.mat','suzuki_yoshida_2011');
return;

% Save data into database
dt = suzuki_yoshida_2011;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));
end

