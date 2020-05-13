function [Idx] = findMinMax(sig,type)

    Idx = [];
    switch type
        case 'min'
            for s = 2:(length(sig)-1)
                if ((sig(s-1)>sig(s)) && (sig(s)<sig(s+1)))
                    Idx = [Idx s];
                end
            end
        case 'max'
            for s = 2:(length(sig)-1)
                if ((sig(s-1)<sig(s)) && (sig(s)>sig(s+1)))
                    Idx = [Idx s];
                end
            end
    end
        
end