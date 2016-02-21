%% Hirasawa~Shimizu, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

hirasawa_shimizu_2013.pmid = 23665193;

phenotypes = {'lactate production'};
treatments = {''};

% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/data_from_S1_PDF.xlsx', 'Mutants');
data.raw(1,:) = [];
data2.genenames = [data.raw(:,1); data.raw(:,4); data.raw(:,7); data.raw(:,10)];
data2.data = [data.raw(:,2:3); data.raw(:,5:6); data.raw(:,8:9); data.raw(:,11:12)];

inds = find(cellfun(@isnumeric, data2.genenames));
data2.genenames(inds) = [];
data2.data(inds,:) = [];

data2.genenames = upper(clean_genename(data2.genenames));

inds = find(~is_genename(data2.genenames) & ~is_orf(data2.genenames));
data2.genenames(inds) = [];
data2.data(inds,:) = [];

[data2.orfs, translated] = translate(data2.genenames);
data2.orfs(~translated) = [];
data2.data(~translated, :) = [];

data2.data = cell2mat(data2.data);
data2.data = nanmean(data2.data,2);

% Load controls
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/data_from_S1_PDF.xlsx', 'CTRL');
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

save('./hirasawa_shimizu_2013.mat','hirasawa_shimizu_2013');

fid = fopen('./hirasawa_shimizu_2013.txt','w');
write_matrix_file(fid, hirasawa_shimizu_2013.orfs, hirasawa_shimizu_2013.ph, hirasawa_shimizu_2013.data);
fclose(fid);

end
