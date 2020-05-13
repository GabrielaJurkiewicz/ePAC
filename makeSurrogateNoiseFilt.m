function [Surrogate] = makeSurrogateNoiseFilt(signal,from,to,fs,Nboot)
    
    Surrogate = zeros(Nboot,length(signal));
    for nb=1:Nboot
        tmp = randn(1,length(signal));
        Surrogate(nb,:) = eegfilt(tmp,fs,from,to); 
    end

end