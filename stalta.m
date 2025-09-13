function [Svalue] = stalta(tseries,ltw,stw,thresh)

% STALTA estimates first arrival of a seismogram based on ratio of
% short-term to long-term averages exceeding some specified 
% threshold value. LTW,STW are window lengths in seconds of 
% long and short term averages, THRESH is threshold and DT is 
% sample interval of time series TSERIES.

    il = ltw; %in samples
    is = stw;
    nt = length(tseries); 
    aseries = abs( hilbert(tseries) );
%     aseries = tseries;
    sra = zeros(nt, 1);
    for ii = il + 1:nt
       lta = mean(aseries(ii - il : ii));
       sta = mean(aseries(ii - is : ii));
       sra(ii) = sta / lta;
       if (lta == 0)
           sra(ii) = 0;
       else
           sra(ii) = sta / lta;
       end
    end
    
%     itm = find(sra > thresh);
%     if ~isempty(itm)
%         itmax = itm(0);
%     end
    
%     indx = itmax;
    Svalue = sra;
    return
    
    