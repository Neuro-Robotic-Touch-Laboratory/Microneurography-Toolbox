tic
indxs = INDEXES;
indxs.dt_mean = nan(size(indxs.t));
for i =1:length(indxs.t)
    indxs.dt_mean(i) = mean(diff(indxs.t_events(indxs.t(i)-window4mean<indxs.t_events & indxs.t_events<indxs.t(i))));
end
toc

tic
INDEXES.dt_mean = arrayfun(windowIntervalMean_fun, INDEXES.t);
toc


isequal (INDEXES.dt_mean,indxs.dt_mean)


