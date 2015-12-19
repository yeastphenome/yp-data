function S = genename2orf(genenames, varargin)

addpath(genpath('~/Laboratory/Datasets/'));
load sgd_141010;

genenames = lower(strtrim(genenames));
[genenames,ia,ic] = unique(genenames);

s = cell(size(genenames));

[~,ind1,ind2] = intersect(genenames, sgd.genenames);
s(ind1) = sgd.orfs(ind2);

% Print out things that couldn't be matched
inds = find(cellfun(@isempty, s));
if ~isempty(inds)
    fprintf('\n\nThese genenames could not be matched:\n')
    for i = 1 : length(inds)
        fprintf('%s\n', genenames{inds(i)});
    end
    fprintf('\n');
end

s(inds) = genenames(inds);

S = s(ic);