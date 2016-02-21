%% Suzuki~Yoshida, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
suzuki_yoshida_2011.pmid = 21601516;

phenotypes = {'glutathione abundance'};
treatments = {'standard'};

% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/hits_data.xlsx');

hits_genes = data.raw(1:end,1);
hits_data = cell2mat(data.raw(1:end,2));

hits_genes = clean_genename(hits_genes);
hits_orfs = translate(hits_genes);

[t,t2] = grpstats(hits_data,hits_orfs,{'gname','mean'});
hits_orfs = t;
hits_data = t2;

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/YKOmatalpha_GSH_list070508.xlsx');
tested_orfs = tested.raw(3:end,3);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(upper(clean_orf(tested_orfs)));

[missing, ix] = setdiff(hits_orfs, tested_orfs);

suzuki_yoshida_2011.orfs = tested_orfs;
suzuki_yoshida_2011.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
suzuki_yoshida_2011.data(ind2,:) = hits_data(ind1,:);

suzuki_yoshida_2011.ph = strcat(phenotypes, '; ', treatments);

save('./suzuki_yoshida_2011.mat','suzuki_yoshida_2011');

fid = fopen('./suzuki_yoshida_2011.txt','w');
write_matrix_file(fid, suzuki_yoshida_2011.orfs, suzuki_yoshida_2011.ph, suzuki_yoshida_2011.data);
fclose(fid);

end
