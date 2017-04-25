%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
kemmeren_holstege_2014.pmid = 24766815;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(kemmeren_holstege_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('read_matrix_file','./raw_data/deleteome_all_mutants_ex_wt_var_controls.txt',3,2);

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data.labels_col(:,1);

% Get the "phenotypes"
hit_ph = data.labels_row(:,2);

% Get the data itself
hit_data = data.data';

% Only keep the M = log2(del/wt)
inds = find(strcmp('M', data.labels_col(:,2)));
hit_strains = hit_strains(inds);
hit_data = hit_data(inds,:);

% Extract the name of the deletion mutants
dels = cellfun(@cell2mat, regexp(hit_strains, '^[a-zA-Z0-9,\(\)]*(?=\-)','match'), 'UniformOutput',0);
inds = find(cellfun(@isempty, dels));
disp(hit_strains(inds));

% A few manual corrections
dels(find(strcmpi('CycC', dels))) = {'SSN8'};
 
hit_genenames = dels;

% Eliminate all white spaces & capitalize
hit_genenames = clean_genename(hit_genenames);

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_genenames);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% If the same gene has more than 2 expression columns, average its values
[hit_ph, hit_data] = grpstats(hit_data', hit_ph, {'gname','mean'});
hit_data = hit_data';

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5658:11769]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
kemmeren_holstege_2014.orfs = hit_strains;
kemmeren_holstege_2014.ph = hit_data_names;
kemmeren_holstege_2014.data = hit_data;
kemmeren_holstege_2014.dataset_ids = hit_data_ids;

%% Save

save('./kemmeren_holstege_2014.mat','kemmeren_holstege_2014');

%% Print out

fid = fopen('./kemmeren_holstege_2014.txt','w');
write_matrix_file(fid, kemmeren_holstege_2014.orfs, kemmeren_holstege_2014.ph, kemmeren_holstege_2014.data);
fclose(fid);

fid = fopen('./data.txt','w');
for i = 1 : length(kemmeren_holstege_2014.ph)
    for j = 1 : length(kemmeren_holstege_2014.orfs)
        if ~isnan(kemmeren_holstege_2014.data(j,i))
            fprintf(fid, '%d\t%s\t%.3f\n', kemmeren_holstege_2014.dataset_ids(i), kemmeren_holstege_2014.orfs{j}, kemmeren_holstege_2014.data(j,i));
        end
    end
end
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(kemmeren_holstege_2014)
end

end

