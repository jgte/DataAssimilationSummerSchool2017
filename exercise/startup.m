addpath(fullfile(getenv('HOME'),'cloud','Common','github','orb','matlab'))
  
dbstop if error
d=dir('*.m');
for i=1:numel(d)
  if strcmp(d(i).name,[mfilename,'.m'])
    continue
  end
  edit(d(i).name)
end
clear d
