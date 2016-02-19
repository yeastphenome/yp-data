%% Fei~Yang, 2008
function FILENAMES = code()
FILENAMES = {};
fei_yang_2008.pmid = 18250201;

phenotypes = {'number of lipid droplets'};
treatments = {'YPD'};


%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('textread','./raw_data/hits.txt','%s %d');

% Get the list of ORFs and the correponding data 
hit_strains = data{1};

% Get the data itself
hit_data = data{2}; % if the dataset is discrete or binary

% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

%% Prepare final dataset

% If the dataset is quantitative:
fei_yang_2008.orfs = hit_strains;
fei_yang_2008.ph = strcat(phenotypes, '; ', treatments);
fei_yang_2008.data = hit_data;

%% Save

save('./fei_yang_2008.mat','fei_yang_2008');

fid = fopen('./fei_yang_2008.txt','w');
write_matrix_file(fid, fei_yang_2008.orfs, fei_yang_2008.ph, fei_yang_2008.data);
fclose(fid);

end

