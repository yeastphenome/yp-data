%% Hirasawa~Shimizu, 2013
% DATA = hirasawa_shimizu_2013
function FILENAMES = code()
FILENAMES = {};

hirasawa_shimizu_2013.pmid = 23665193;

phenotypes = {'lactate production'};
treatments = {''};

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/data_from_S1_PDF.xlsx', 'Mutants');
data.raw(1,:) = [];
data2.genenames = [data.raw(:,1); data.raw(:,4); data.raw(:,7); data.raw(:,10)];
data2.data = [data.raw(:,2:3); data.raw(:,5:6); data.raw(:,8:9); data.raw(:,11:12)];

inds = find(cellfun(@isnumeric, data2.genenames));
data2.genenames(inds) = [];
data2.data(inds,:) = [];

data2.genenames = cellfun(@strtrim, data2.genenames,'UniformOutput',0);

data2.orfs = genename2orf(data2.genenames,'noannot');

% Adjustments
data2.orfs(strcmp('ARR4', data2.orfs)) = {'YDL100C'};
data2.orfs(strcmp('FMP12', data2.orfs)) = {'YHL021C'};
data2.orfs(strcmp('FMP13', data2.orfs)) = {'YKR016W'};
data2.orfs(strcmp('FMP14', data2.orfs)) = {'YPL099C'};
data2.orfs(strcmp('FMP22', data2.orfs)) = {'YHR198C'};
data2.orfs(strcmp('FMP24', data2.orfs)) = {'YMR115W'};
data2.orfs(strcmp('FMP26', data2.orfs)) = {'YJR080C'};
data2.orfs(strcmp('FMP29', data2.orfs)) = {'YER080W'};
data2.orfs(strcmp('FMP34', data2.orfs)) = {'YHR199C'};
data2.orfs(strcmp('FMP36', data2.orfs)) = {'YDR493W'};
data2.orfs(strcmp('FMP38', data2.orfs)) = {'YOR205C'};
data2.orfs(strcmp('FMP39', data2.orfs)) = {'YMR157C'};
data2.orfs(strcmp('FMP50', data2.orfs)) = {'YKR027W'};
data2.orfs(strcmp('FMP51', data2.orfs)) = {'YBR262C'};
data2.orfs(strcmp('FUN34', data2.orfs)) = {'YNR002C'};
data2.orfs(strcmp('GRD19', data2.orfs)) = {'YOR357C'};
data2.orfs(strcmp('MDM39', data2.orfs)) = {'YGL020C'};
data2.orfs(strcmp('MSU1', data2.orfs)) = {'YMR287C'};
data2.orfs(strcmp('PRP12', data2.orfs)) = {'YMR302C'};
data2.orfs(strcmp('RCS1', data2.orfs)) = {'YGL071W'};
data2.orfs(strcmp('RMD7', data2.orfs)) = {'YER083C'};
data2.orfs(strcmp('SOY1', data2.orfs)) = {'YBR194W'};
data2.orfs(strcmp('ZMS1', data2.orfs)) = {'YJR127C'};

data2.orfs = upper(data2.orfs);

inds = find(~strncmp('Y', data2.orfs,1));
data2.orfs(inds) = [];
data2.genenames(inds) = [];
data2.data(inds,:) = [];

data2.data = cell2mat(data2.data);
data2.data = nanmean(data2.data,2);

% Load controls
[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/data_from_S1_PDF.xlsx', 'CTRL');
data.raw(1,:) = [];
ctrl_data = [data.raw(:,2:3); data.raw(:,5:6); data.raw(:,8:9); data.raw(:,11:12)];
ctrl_data = reshape(ctrl_data,[],1);
inds = find(~cellfun(@isnumeric, ctrl_data));
ctrl_data(inds) = [];
ctrl_data = cell2mat(ctrl_data);
ctrl_data = nanmean(ctrl_data);

data2.data_norm = data2.data ./ ctrl_data;

% Average replicates
[t,t2] = grpstats(data2.data_norm, data2.orfs, {'gname','mean'});
hirasawa_shimizu_2013.orfs = t;
hirasawa_shimizu_2013.data = t2;
hirasawa_shimizu_2013.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'hirasawa_shimizu_2013.mat'],'hirasawa_shimizu_2013');
return;

% Save data into database
dt = hirasawa_shimizu_2013;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));



end

