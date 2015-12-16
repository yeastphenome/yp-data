%% Hancock~Lopes, 2006
function FILENAMES = code()
FILENAMES = {};
hancock_lopes_2006.pmid = 16582425;

phenotypes = {'Opi-'};
treatments = {''};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/mat_alpha_061101.xlsx', 'mat_alpha_061101');
tested_orfs = tested.raw(4:end,2);

inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];

tested_orfs = regexprep(tested_orfs, '\W','');

inds = find(cellfun(@isempty, regexpi(tested_orfs,'Y[A-P][RL][0-9]{3}[CW](-[ABC])*')));
tested_orfs(inds) = [];

tested_orfs = unique(upper(tested_orfs));

% Load data
fid = fopen('raw_data/hits_genenames.txt');
hits_genenames = textscan(fid, '%s');
hits_genenames = hits_genenames{1};
fclose(fid);

hits_orfs = genename2orf(hits_genenames,'noannot');

hits_orfs(strcmp(hits_orfs,'vma12')) = {'YKL119C'};
hits_orfs(strcmp(hits_orfs,'ada3')) = {'YDR176W'};
hits_orfs(strcmp(hits_orfs,'rcs1')) = {'YGL071W'};
hits_orfs(strcmp(hits_orfs,'cwh36')) = {'YCL005W-A'};
hits_orfs(strcmp(hits_orfs,'rmd7')) = {'YER083C'};

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];

hits_orfs = regexprep(hits_orfs,'\W','');

inds = find(cellfun(@isempty, regexpi(hits_orfs,'Y[A-P][RL][0-9]{3}[CW](-[ABC])*')));
hits_orfs(inds) = [];

hits_orfs = unique(upper(hits_orfs));

hits_scores = ones(length(hits_orfs),1);

hancock_lopes_2006.orfs = tested_orfs;
hancock_lopes_2006.data = zeros(length(tested_orfs), length(phenotypes));
[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
hancock_lopes_2006.data(ind2,1) = hits_scores(ind1,1);

hancock_lopes_2006.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'hancock_lopes_2006.mat'],'hancock_lopes_2006');
return;

% Save data into database
dt = hancock_lopes_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


end

