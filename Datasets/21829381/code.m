%% Fei~Yang, 2011
function FILENAMES = code()
FILENAMES = {};
fei_yang_2011.pmid = 21829381;

phenotypes = {'cells with superlarge lipid droplets'};
treatments = {'SC'; 'YPD'};


%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('read_matrix_file','./raw_data/hits.txt');

% Get the list of ORFs and the correponding data 
hit_strains = data.labels_row(:,1);

% Get the data itself
hit_data = data.data; % if the dataset is discrete or binary

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% Normalize to WT
indWt = find(strcmp('WT', hit_strains));
hit_data = hit_data ./ repmat(hit_data(indWt,:),length(hit_strains),1);

hit_strains(indWt) = [];
hit_data(indWt,:) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});


%% Prepare final dataset

fei_yang_2011.orfs = hit_strains;
fei_yang_2011.ph = strcat(phenotypes, '; ', treatments);
fei_yang_2011.data = hit_data;

%% Save

save('./fei_yang_2011.mat','fei_yang_2011');

fid = fopen('./fei_yang_2011.txt','w');
write_matrix_file(fid, fei_yang_2011.orfs, fei_yang_2011.ph, fei_yang_2011.data);
fclose(fid);

end

