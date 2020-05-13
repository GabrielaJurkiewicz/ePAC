function [spectrumInterp] = fitBackgroundSpectrum(spectrum,f)

    [IdxMin] = findMinMax(spectrum,'min');
    IdxMin = [1 IdxMin length(f)];
    IdxMin = sort(unique(IdxMin));
    spectrumInterp = pchip(f(IdxMin),spectrum(IdxMin),f);

    [p1] = findMinMax(diff(spectrum),'min');
    [p2] = findMinMax(diff(spectrum),'max');
    PP = sort([p1,p2]);
    for ff = 1:length(f)
        if spectrumInterp(ff)>spectrum(ff)
            newP = min(PP(PP>ff));
            if ~isempty(newP)
                if sum(IdxMin==newP)==0
                    spInterp = pchip(f(sort([IdxMin newP])),spectrum(sort([IdxMin newP])),f);
                    s1 = sum(spectrumInterp(spectrumInterp>spectrum)-spectrum(spectrumInterp>spectrum));
                    s2 = sum(spInterp(spInterp>spectrum)-spectrum(spInterp>spectrum));
                    if s2<s1
                        IdxMin = sort([IdxMin newP]);
                        spectrumInterp = spInterp;
                    end
                end
            end
        end
    end

end