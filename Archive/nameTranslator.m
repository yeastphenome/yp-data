% This function will take an array of names in any format (DBID, standard,
% systematic, or by name) and will return an array in the format that you 
% want with all possible changes (DBID, standard, systematic, or by name)

function newNames = nameTranslator(typehave, typeneed, oldNames)

load('toSystematic.mat');

switch upper(typehave)
    case 'DBID'
        inds = cellfun(@(x) find(ismember(toSystematic.primaryDBID, x)), oldNames, 'UniformOutput', false);
        newNames = retrieveName();
    case 'SYSTEMATIC'
        inds = cellfun(@(x) find(ismember(toSystematic.systematicName, x)), oldNames, 'UniformOutput', false);
        newNames = retrieveName();
    case 'STANDARD'
        inds = cellfun(@(x) find(ismember(toSystematic.standardName, x)), oldNames, 'UniformOutput', false);
        newNames = retrieveName();
    case 'NAME'
        inds = cellfun(@(x) find(ismember(toSystematic.name, x)), oldNames, 'UniformOutput', false);
        newNames = retrieveName();
    otherwise
        error('Typehave must be DBID, systematic, standard, or name');
end
    function newNames = retrieveName()
        indx = find(cellfun(@isempty, inds));
        indx2 = find(~cellfun(@isempty, inds));
        inds(indx) = [];
        newNames = oldNames;
        switch upper(typeneed)
            case 'DBID'
                newNames(indx2) = cellfun(@(x) toSystematic.primaryDBID{x}, inds, 'UniformOutput', false);
            case 'SYSTEMATIC'
                newNames(indx2) = cellfun(@(x) toSystematic.systematicName{x}, inds, 'UniformOutput', false);
            case 'STANDARD'
                newNames(indx2) = cellfun(@(x) toSystematic.standardName{x}, inds, 'UniformOutput', false);
            case 'NAME'
                newNames(indx2) = cellfun(@(x) toSystematic.name{x}, inds, 'UniformOutput', false);
            otherwise
                error('Typeneed must be DBID, systematic, standard, or name');
        end
end

return;
end