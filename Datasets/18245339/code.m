%% Abe~Minegishi, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
abe_minegishi_2008.pmid = 18245339;

phenotypes = {'growth (culture turbidity)'};
treatments = {'high pressure';'low temperature'};

% Load tested strains
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/mat_alpha_041902.xlsx', 'mat_alpha_041902.txt');

tested_orfs = tested.raw(4:end,3);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(upper(clean_orf(tested_orfs)));

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Table1_Abe_Genetics.xlsx', 'Table 1 (2)');

hits_orfs = data.raw(7:end,3);
hits_data = data.raw(7:end,[26 31]);

inds = find(cellfun(@isempty, hits_orfs) | cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_orfs = clean_orf(hits_orfs);

hits_data = cell2mat(hits_data);
hits_data = hits_data/100;  % transform percent into fractions

[missing, ix] = setdiff(hits_orfs, tested_orfs);    % nothing missing

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(hits_data, hits_orfs, {'gname','mean'});

abe_minegishi_2008.orfs = tested_orfs;
abe_minegishi_2008.data = zeros(length(tested_orfs),2);

[~,ind1,ind2] = intersect(t, tested_orfs);
abe_minegishi_2008.data(ind2,:) = t2(ind1,:);

abe_minegishi_2008.ph = strcat(phenotypes, '; ', treatments);

save('./abe_minegishi_2008.mat','abe_minegishi_2008');

fid = fopen('./abe_minegishi_2008.txt','w');
write_matrix_file(fid, abe_minegishi_2008.orfs, abe_minegishi_2008.ph, abe_minegishi_2008.data);
fclose(fid);

end
