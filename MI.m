function [mi] = MI(Phase,Amp,phaseBins)

    nbin = length(phaseBins);  
    MeanAmp = zeros(1,nbin); 
    for j = 1:nbin   
        I = ((Phase <  (phaseBins(j)+(2*pi/nbin))) & (Phase >=  phaseBins(j)));
        MeanAmp(j) = mean(Amp(I)); 
    end
    
    % quantifying the amount of amp modulation by means of a normalized entropy index (Tort et al PNAS 2008):
    mi = (log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);

end