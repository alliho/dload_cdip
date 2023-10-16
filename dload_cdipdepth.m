function [depth] = dload_cdipdepth(cdipid)
% -------------------------------------------------------------------------
% dload_cdipdepth  Gets depth of CDIP buoy
% -------------------------------------------------------------------------
%   Syntax:
%      [depth] = dload_cdipdepth(cdipid)
%
%   Inputs:
%      CDIPID 
% 
%   Output:
%      [depth]
% 
% Updated as of 01-23-2020 by Alli Ho
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% function set up 
idstr = cdipid;

if isnumeric(idstr)
    idstr = num2str(idstr);
    if length(idstr)<3
        idstr = ['0' idstr];
    end
end


%% load
try
    url = ['http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/' idstr 'p1_rt.nc.ascii?metaWaterDepth'];
    dat = urlread(url); dat = strsplit(dat,'\n');
    idx = find(contains(dat, '---------------------------------------------')); 
    tmp = strsplit(dat{idx+1}, ',');
    depth = str2num(tmp{2});

catch % historic only
    url = ['http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/' idstr 'p1/' idstr 'p1_historic.nc.ascii?metaWaterDepth'];
    dat = urlread(url); dat = strsplit(dat,'\n');
    idx = find(contains(dat, '---------------------------------------------')); 
    tmp = strsplit(dat{idx+2}, ',');
    depths = str2num([tmp{1:end}]);
%     depth = nanmean(depths);
    depth = depths;
%     https://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/240p1/240p1_historic.nc.ascii?metaWaterDepth%5B0:1:5%5D

end


end
