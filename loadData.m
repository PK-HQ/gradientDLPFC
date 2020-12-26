function [datasets, bins_overlap, regional, sessions, trials, m, st] = loadData(datapath, timeBins)
load(datapath)
datasets=datasets(:,:,timeBins);
bins_overlap=bins_overlap(:,timeBins);
end