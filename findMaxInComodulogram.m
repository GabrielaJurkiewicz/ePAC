function [fmx,fwd,fstart,fstop,fspmx,spmx] = findMaxInComodulogram(comod,SA,AS,f,fA,w,Fs)

    fmx = [];
    fwd = [];
    fstart = [];
    fstop =[];
    fspmx = [];
    spmx = [];
    if sum(comod>0)>0
        ind = (comod>0)';
        ind = ind(2:end)-ind(1:end-1);
        start = find(ind==1);
        stop = find(ind==-1);
        if isempty(start)
            start = [1];
        else
            start = start+1;
        end
        if isempty(stop)
            stop = [length(comod)];
        else
            stop = stop;
        end
        if stop(1)<start(1)
            start = [1 start]; 
        end
        if start(end)>stop(end)
            stop = [stop length(comod)];
        end
        for ss=1:length(start)
            [~,idMaxPAC] = max(comod(start(ss):stop(ss)));
            if ~isempty(idMaxPAC)
                fmx = [fmx fA(idMaxPAC+start(ss)-1)];
                wdt = findWaveletFreqWidth(fA(idMaxPAC+start(ss)-1),w,Fs);
                fwd = [fwd wdt];
                fstart = [fstart fA(start(ss))];
                fstop = [fstop fA(stop(ss))];
                
                idf = (f>=fA(start(ss)))&(f<=fA(stop(ss)));
                idfMAX = ((f>=(fA(idMaxPAC+start(ss)-1)-wdt/2)) & (f<=(fA(idMaxPAC+start(ss)-1)+wdt/2)));
                ID = find(idfMAX|idf);
                as = AS(ID);
                sa = SA(ID);
                tmpas = sum(as(as>sa));
                sa = sum(sa(sa>as));
                as = tmpas;
                
                if (sa>=as || sa>(0.999999*as))
                    [~,idsa] = findpeaks(SA(ID));
                    [~,idM] = max(SA(ID));
                    if sum(idsa==idM)>0
                        idsa = idsa + ID(1) - 1;
                        [~,idmsa] = max(SA(idsa));
                        idsa = idsa(idmsa);
                        fspmx = [fspmx f(idsa)];
                        spmx = [spmx SA(idsa)];
                    else
                        fspmx = [fspmx f(idM + ID(1) - 1)];
                        spmx = [spmx SA(idM + ID(1) - 1)];
                    end
                else
                    [~,idas] = findpeaks(AS(ID));
                    [~,idM] = max(AS(ID));
                    if sum(idas==idM)>0
                        idas = idas + ID(1) - 1;
                        [~,idmas] = max(AS(idas));
                        idas = idas(idmas);
                        fspmx = [fspmx f(idas)];
                        spmx = [spmx AS(idas)];
                    else
                        fspmx = [fspmx f(idM + ID(1) - 1)];
                        spmx = [spmx AS(idM + ID(1) - 1)];
                    end
                end
            end
        end
    end

end