%% Nyswaner~Garfinkel,2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
nyswaner_garfinkel_2008.pmid = 18202368;

phenotypes = {'number of cells (His+) to monitor the frequency of Ty1 transposition'};
treatments = {''};

%% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Matalphakos counted.xlsx', 'Sheet2');
tested_orfs = tested.raw(6:end,2);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = clean_orf(tested_orfs);

% Fix typos
tested_orfs(find(strcmp('YLR228', tested_orfs))) = {'YLR228C'};
tested_orfs(find(strcmp('YMR062', tested_orfs))) = {'YMR062C'};
tested_orfs(find(strcmp('YYKL138C', tested_orfs))) = {'YKL138C'};

inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));

tested_orfs(inds) = [];


%% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/genetics.107.082602-9.xlsx', 'Sheet1');
hits_orfs = data.raw(:,1);

inds = find(cellfun(@isempty, hits_orfs) | cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];

hits_orfs = clean_orf(hits_orfs);

inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));

hits_orfs(inds) = [];

missing = setdiff(hits_orfs, tested_orfs);

% Adding 3 ORFs to the list of tested
tested_orfs = [tested_orfs; missing];

nyswaner_garfinkel_2008.orfs = tested_orfs;
nyswaner_garfinkel_2008.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(tested_orfs, hits_orfs);
nyswaner_garfinkel_2008.data(ind1,1) = 1;

nyswaner_garfinkel_2008.ph = strcat(phenotypes, '; ', treatments);

save('./nyswaner_garfinkel_2008.mat','nyswaner_garfinkel_2008');

fid = fopen('./nyswaner_garfinkel_2008.txt','w');
write_matrix_file(fid, nyswaner_garfinkel_2008.orfs, nyswaner_garfinkel_2008.ph, nyswaner_garfinkel_2008.data);
fclose(fid);

end

