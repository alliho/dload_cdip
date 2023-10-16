function [name] = dload_cdipname(cdipid)
% -------------------------------------------------------------------------
% dload_cdipname  Gets formal name of CDIP buoy
% -------------------------------------------------------------------------
%   Syntax:
%      [name] = dload_cdipname(cdipid)
%
%   Inputs:
%      CDIPID 
% 
%   Output:
%      [name]
% 
% Updated as of 02-13-2023 by Alli Ho
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



urlstr = ['https://cdip.ucsd.edu/m/products/?stn=' idstr];
rawdata = urlread(urlstr); data = strsplit(rawdata,'\n');
dstart = find(contains(data,[idstr ' - ']));
var = strsplit(data{dstart}, '-'); var = var{2}; var = var(2:end-1);
name = var;



end


