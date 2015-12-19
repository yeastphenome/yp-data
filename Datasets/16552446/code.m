%% Gatbonton~Bedalov, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

gatbonton_bedalov_2006.source = {'main PDF'};
gatbonton_bedalov_2006.downloaddate = {'2014-03-10'};
gatbonton_bedalov_2006.pmid = 16552446;

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/gatbonton_bedalov_2006_hits.xlsx', 'Sheet1');

[hits_orfs, translated] = translate(data.raw(:,1));
hits_orfs(~translated) = [];
data.raw(~translated,:) = [];

hits_scores = cell2mat(data.raw(:,3));

% Converting the scores such that:
% a) 1 = weakest phenotype, 3 = strongest phenotype
% b) long telomeres = positive scores, short telomeres = negative scores
hits_scores = abs(hits_scores-4);
inds = find(strcmp('S', data.raw(:,2)));
hits_scores(inds) = -hits_scores(inds);

phenotypes = {'telomere length'};
treatments = {''};


% Load tested genes
[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/genelist_altered_020806.xlsx', 'mat alpha copy.txt');

tested_orfs = data.raw(2:end,1);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs(strcmp('YYKL138C', tested_orfs)) = {'YKL138C'};
tested_orfs = unique(upper(cleanOrf(tested_orfs)));

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
gatbonton_bedalov_2006.orfs = tested_orfs;
gatbonton_bedalov_2006.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(gatbonton_bedalov_2006.orfs, hits_orfs);
gatbonton_bedalov_2006.data(ind1,:) = hits_scores(ind2,:);

gatbonton_bedalov_2006.ph = [strcat(phenotypes{1}, '; ', treatments)];


a = mfilename('fullpath');
a = a(1:end-4);
save([a,'gatbonton_bedalov_2006.mat'],'gatbonton_bedalov_2006');
return;

% Save data into database
dt = gatbonton_bedalov_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


end

