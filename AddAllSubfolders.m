% AddAllSubfolders - dodaje wszystkie podfoldery do path MATLAB
% Wyklucza katalogi systemowe (.git) i ukryte katalogi konfiguracyjne
allpaths = genpath(pwd);
paths = strsplit(allpaths, pathsep);
paths = paths(~cellfun(@(p) contains(p, '.git') || startsWith(p, '.'), paths));
addpath(strjoin(paths, pathsep));
