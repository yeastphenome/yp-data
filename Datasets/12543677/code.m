%% Blackburn~Avery, 2003
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
% TESTED = not available

blackburn_avery_2003.source = {'manuscript PDF'};
blackburn_avery_2003.downloaddate = {'2014-02-12'};
blackburn_avery_2003.pmid = 12543677;

[FILENAMES{end+1}, data.raw] = readdata('xlsread','./raw_data/blackburn_avery_2003_data.xlsx', 'data.txt');

phenotypes = {'growth (MIC)'};
treatments = data.raw(1,2:8)';
data.raw(1,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,1) = cleanOrf(data.raw(:,1));

% Eliminate everything that doesn't look like an ORF
inds = find(~isorf(data.raw(:,1)));
data.raw(inds,:) = [];

% Replace 'Inf' with Inf
t = data.raw(:,2:end);
t(~cellfun(@isnumeric, t)) = {Inf};

blackburn_avery_2003.orfs = upper(data.raw(:,1));
blackburn_avery_2003.data = cell2mat(t);

blackburn_avery_2003.ph = [strcat(phenotypes{1}, '; ', treatments)];

save('./blackburn_avery_2003.mat','blackburn_avery_2003');
return;

% Save data into database
dt = blackburn_avery_2003;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


end

