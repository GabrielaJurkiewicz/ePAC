function [wd] = findWaveletFreqWidth(f,w,Fs)

    T = 1;
    s = T*f/(2*w);
    psi = fft(tf_morlet(Fs,w,s,true));
    psi = psi./sqrt(sum(abs(psi).^2));
    psi = sqrt(psi.*conj(psi));
    psi = 2*psi(1:floor(Fs/2)+1);
    [~,idm,wd,~] = findpeaks(psi);
    freq = Fs/2*linspace(0,1,Fs/2+1);
    wd = wd((freq(idm)>=(f-1))&(freq(idm)<=(f+1)));
%     plot(freq,psi)

end