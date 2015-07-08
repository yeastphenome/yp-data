%% Ohya~Morishita, 2005
% DATA = ohya_morishita_2005
function FILENAMES = code()
FILENAMES = {};


ohya_morishita_2005.source = {'http://scmd.gi.k.u-tokyo.ac.jp/datamine/pnas/mutant_analysis_2011_10_20.tab'};
ohya_morishita_2005.downloaddate = {'2013-03-12'};
ohya_morishita_2005.pmid = 16365294;

phenotypes = {'Growth, flow cytometry'};
treatments = {''};

[FILENAMES{end+1}, param.raw] = dataread('xlsread','Datasets/Phenotypes/2005_Ohya~Morishita/parameter_names.xlsx', 'Sheet1');

[FILENAMES{end+1}, data] = dataread('importdata','Datasets/Phenotypes/2005_Ohya~Morishita/mutant_analysis_2011_10_20.tab');
data.orfs = data.textdata(2:end,1);
data.params = data.textdata(1,2:end)';

[t,ind1,ind2] = intersect(data.params, param.raw(:,1));
data.params_name(ind1,1) = param.raw(ind2,2);


% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.orfs));
data.orfs(inds) = [];
data.data(inds,:) = [];

inds = setdiff(1:length(data.orfs),strmatch('Y', data.orfs));
data.orfs(inds) = [];
data.data(inds,:) = [];

% Eliminate white spaces before/after ORF
data.orfs = cellfun(@strtrim, data.orfs,'UniformOutput',0);
data.orfs = upper(data.orfs);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data.data, data.orfs, {'gname','mean'});
ohya_morishita_2005.orfs = t;
ohya_morishita_2005.data = t2;
ohya_morishita_2005.ph = data.params_name;

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'ohya_morishita_2005.mat'],'ohya_morishita_2005');
return;


% NEED TO SEND EMAIL TO YOSHI OHYA WITH QUESTIONS
% % Save data into database
% dt = ohya_morishita_2005;
% 
% datasets = get_datasets_for_paper(dt);
% datasets_ids = zeros(length(datasets),1);
% datasets_names = cell(length(datasets),3);
% for i = 1 : length(datasets)
%     datasets_ids(i,1) = datasets(i).id;
%     datasets_names{i,1} = datasets(i).name;
%     if isempty(datasets(i).reporter)
%         datasets_names{i,2} = '';
%     else
%         datasets_names{i,2} = datasets(i).reporter;
%     end
%     datasets_names{i,3} = datasets(i).conditionset;
% end
% 
% [~,database_ix] = sortrows(datasets_names,[1 2 3]);
% [~,ph_ix] = sort(dt.ph);
% 
% % Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
% datasets_names(database_ix,:)
% dt.ph(ph_ix)
% 
% insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

end

