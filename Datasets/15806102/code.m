%% Giorgini~Muchowski, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

giorgini_muchowski_2005.pmid = 15806102;

[FILENAMES{end+1}, hits] = read_data('textread','./raw_data/giorgini_muchowski_2005_hits.txt', '%s');

phenotypes = {'growth (pooled CFU)'};
treatments = {'Htt103Q'};

% Eliminate white spaces before/after gene names
hits(:,1) = clean_genename(hits);

% Translate genenames to ORF
hits_orfs = translate(hits);

% Hits = LOF suppressors of Htt103Q sensitivity, so positive score.
scores = ones(length(hits_orfs),1);

% Load tested genes
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat_a_obs_v2(1).0.xlsx', 'DATA');
inds = find(cellfun(@isnumeric, tested.raw(:,2)));
tested.raw(inds,:) = [];

tested_orfs = unique(upper(tested.raw(2:end,2)));
tested_orfs(strcmp(tested_orfs,'YLR287-A')) = {'YLR287C-A'};

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);
tested_orfs = [tested_orfs; missing];

% Create dataset
giorgini_muchowski_2005.orfs = tested_orfs;
giorgini_muchowski_2005.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(giorgini_muchowski_2005.orfs, hits_orfs);
giorgini_muchowski_2005.data(ind1,:) = scores(ind2,:);

giorgini_muchowski_2005.ph = [strcat(phenotypes{1}, '; ', treatments)];


save('./giorgini_muchowski_2005.mat','giorgini_muchowski_2005');
return;

% Save data into database
dt = giorgini_muchowski_2005;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

fid = fopen('./giorgini_muchowski_2005.txt','w');
write_matrix_file(fid, giorgini_muchowski_2005.orfs, giorgini_muchowski_2005.ph, giorgini_muchowski_2005.data);
fclose(fid);

end

