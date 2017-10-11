function [selectedindex,corrmat]=xspecmatStats(freqarraystack,freqrange)
% INPUT:
% freqarraystack: complex cross-spectra matrix of one or two stations; 
%                 nwindows*nfreqstore matrix;
% freqrange:      1*nfreqstore vector; corresponds to all discrete frequencies stored in freqarraystack matrix
%
% RETURN: 
% selectedindex: boolean vector of selected rows in freqarraystack where 0==outliers; 1==inliers
% corrmat:       correlation matrix of neighboring frequency power spectral values

% Liu, X., Ben-Zion, Y., & Zigone, D. (2016). 
% Frequency domain analysis of errors in cross-correlations of ambient seismic noise. Geophysical Journal International, 207(3), 1630-1652.

% Liu, X., & Ben-Zion, Y. (2016). Estimating correlations of neighbouring frequencies 
% in ambient seismic noise. Geophysical Journal International, 206(2), 1065?1075. doi:10.1093/gji/ggw196

% by Xin Liu
% liuxin@stanford.edu
% USC, Dec 2015
% updated at Stanford, Oct 2016
% Non-commercial use only. Commercial use is prohibited without
% permission from the author
% All rights reserved.

% DF=0.05;% Hz
% freqrange=0.05:DF:0.6-DF;% Hz

% No. of windows
nwindows=size(freqarraystack,1);

DF=freqrange(2)-freqrange(1);

freqlist=( freqrange-freqrange(1) )/DF+1;% index of array
nfreqstore=length(freqlist);

boolmat=false(nwindows,nfreqstore);

% outlier threshold:
threshold=3;% 3 times MAD


% noutliers=round(niter*percentoutlier);% integer
for ifreqdistri=1:nfreqstore
    ifreqdistri
%     ifreqdistri=39;% used in the draft fig
    %freqsample=outliers((freqarraystack(:,ifreqdistri)),100);
    indbool=madoutlier(freqarraystack(:,ifreqdistri),threshold);
    freqallsamp=freqarraystack(:,ifreqdistri);
    selectcol=~indbool;
    boolmat(:,ifreqdistri)=selectcol;
    nsamps=sum(selectcol);
    % mean of selected windows:
    mx=mean(freqallsamp(selectcol));
    %%% COMPUTE 2ND MODMENT VARS
    % standard error:
    stderr=std(real(freqallsamp(selectcol)))/sqrt(nsamps);
    

    
end


%% SORT OUT THE SAME OUTLIERS ACROSS FREQ
percentoutlier = 0.08;
goodcells=sum(boolmat,2);

selectedindex=( (nfreqstore-goodcells)<nfreqstore*percentoutlier);% 0.08 is like the tolerance percentage for outliers;
nselected=sum(selectedindex);


corrmat=corrmatesti(freqarraystack(selectedindex,:));


end

function corrmat=corrmatesti(freqarraystackselect)
% THE CORR MATRIX
% Liu, X., Ben-Zion, Y., & Zigone, D. (2016). 
% Frequency domain analysis of errors in cross-correlations of ambient seismic noise. Geophysical Journal International, 207(3), 1630-1652.

%%% by Xin Liu
%%% USC, Oct 2015

corrmat=(corr(real(freqarraystackselect )));% WORKS

end

function indbool=madoutlier(colvec,threshold)
% INPUT: colvec is REAL column vector
% RETURN: indbool where 1==outliers
% %%% outlier removal algorithm based on MAD
% Liu, X., Ben-Zion, Y., & Zigone, D. (2016). 
% Frequency domain analysis of errors in cross-correlations of ambient seismic noise. Geophysical Journal International, 207(3), 1630-1652.

%%% by Xin Liu
%%% USC, Oct 2015

% threshold=3;% 3 times MAD

medRe=median(colvec);

madcolvec=mad(colvec);

indbool=abs(colvec-medRe)./madcolvec>threshold;

end
