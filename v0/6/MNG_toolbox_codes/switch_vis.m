function switch_vis(filename)
tic
fgr =load(filename,'-mat' ,'hgM_070000' );
toc
tic
fgr.hgM_070000.GraphicsObjects.Format3Data.CreateFcn  = 'set(gcf,''visible'',''on'')';
toc
tic
save(filename,'-struct','fgr','-append')
toc
end

