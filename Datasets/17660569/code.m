%% Wilson~van Hoof, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
wilson_vanhoof_2007.pmid = 17660569;

phenotypes = {'growth (suppression of mRNA decay)'};
treatments = {'SC -His'};


%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/hit_list.xlsx', 'Sheet1');

% Get the list of ORFs
hit_strains = data(3:end,1);

% Get the data itself
hit_data = data(3:end,3);
    
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds) = [];

hit_data = cellfun(@length, hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

%% Load tested strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/deletion database.xlsx', 'Sheet1');

% Get the list of ORFs
tested_strains = data(2:end,2);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% Fix typos
tested_strains(strcmp('YPL006', tested_strains)) = {'YPL006W'};
tested_strains(strcmp('YDR007', tested_strains)) = {'YDR007W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds)); 

tested_strains(inds) = [];

tested_strains = unique(tested_strains);

[missing,~] = setdiff(hit_strains, tested_strains); % 0 found.


%% Prepare final dataset

wilson_vanhoof_2007.orfs = tested_strains;
wilson_vanhoof_2007.ph = strcat(phenotypes, '; ', treatments);
wilson_vanhoof_2007.data = zeros(length(wilson_vanhoof_2007.orfs), length(wilson_vanhoof_2007.ph));

[~,ind1,ind2] = intersect(hit_strains, tested_strains);
wilson_vanhoof_2007.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./wilson_vanhoof_2007.mat','wilson_vanhoof_2007');

%% Print out

fid = fopen('./wilson_vanhoof_2007.txt','w');
write_matrix_file(fid, wilson_vanhoof_2007.orfs, wilson_vanhoof_2007.ph, wilson_vanhoof_2007.data);
fclose(fid);

end

