%% Hellauer~Turcotte, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

hellauer_turcotte_2005.source = {'main PDF'};
hellauer_turcotte_2005.downloaddate = {'2014-03-04'};
hellauer_turcotte_2005.pmid = 16061773;

[FILENAMES{end+1}, DATA] = read_data('textread','./raw_data/hits_all.txt', '%s %d');

hits_orfs = DATA{1};
scores = DATA{2};

phenotypes = {'growth (colony size)'};
treatments = {'tirapazamine'};

% Eliminate white spaces before/after ORF
hits_orfs(:,1) = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

% Load tested genes
[FILENAMES{end+1}, tested_orfs] = read_data('textread','./raw_data/ORF.txt', '%s');

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
hellauer_turcotte_2005.orfs = tested_orfs;
hellauer_turcotte_2005.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(hellauer_turcotte_2005.orfs, hits_orfs);
hellauer_turcotte_2005.data(ind1,:) = scores(ind2,:);

hellauer_turcotte_2005.ph = [strcat(phenotypes{1}, '; ', treatments)];


save('./hellauer_turcotte_2005.mat','hellauer_turcotte_2005');
return;

% Save data into database
dt = hellauer_turcotte_2005;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

fid = fopen('./hellauer_turcotte_2005.txt','w');
write_matrix_file(fid, hellauer_turcotte_2005.orfs, hellauer_turcotte_2005.ph, hellauer_turcotte_2005.data);
fclose(fid);

end

