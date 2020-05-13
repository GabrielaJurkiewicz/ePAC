function [X,Y,Colors,L_rel,L_ambig] = findBoundaries(fP,fA,comodulogram,CouplingOrigins)

    new_fP = fP(1):(fP(2)-fP(1))/4:fP(end);
    new_fA = fA(1):(fA(2)-fA(1))/4:fA(end);
    mask = (comodulogram~=0);
    new_mask_rel = zeros(length(new_fA),length(new_fP));
    new_mask_ambig = zeros(length(new_fA),length(new_fP));
    for fp = 1:length(fP)
        for fa = 1:length(fA)
            if mask(fa,fp)==1
                ifa = find(new_fA==fA(fa));
                ifp = find(new_fP==fP(fp));
                sta = ifa-2;
                spa = ifa+2;
                stp = ifp-2;
                spp = ifp+2;
                if fa==length(fA)
                    spa = ifa;
                elseif fa==1
                    sta = ifa;
                end
                if fp==length(fP)
                    spp = ifp;
                elseif fp==1
                    stp = ifp;
                end
                if CouplingOrigins(fa,fp)==1
                    new_mask_rel(sta:spa,stp:spp) = 1;
                else
                    new_mask_ambig(sta:spa,stp:spp) = -1;
                end
            end
        end
    end
    [B_rel,~] = bwboundaries(new_mask_rel,'noholes');
    [B_ambig,~] = bwboundaries(new_mask_ambig,'noholes');
    [~,L_rel] = bwboundaries(comodulogram.*CouplingOrigins,'noholes');
    [~,L_ambig] = bwboundaries(comodulogram-comodulogram.*CouplingOrigins,'noholes');
    base_rel = [[24 150 120]/255; [84,204,182]/255; [169,247,233]/255];
    base_ambig = [[114,77,142]/255; [155,115,188]/255; [221,199,237]/255];
    if ~isempty(B_rel)
        if length(B_rel)>3
            old = linspace(0, 1, length(base_rel));
            new = linspace(0, 1, length(B_rel));
            colors_rel = zeros(length(B_rel), 3);
            for i=1:3
                colors_rel(:,i) = min(max(interp1(old, base_rel(:,i), new)', 0), 1);
            end
        elseif length(B_rel)==3
            colors_rel = base_rel(:,:);
        elseif length(B_rel)==2
            colors_rel = base_rel([1,3],:);
        else
            colors_rel = base_rel(1,:);
        end
    else
        colors_rel = [];
    end
    if ~isempty(B_ambig)
        if length(B_ambig)>3
            old = linspace(0, 1, length(base_ambig));
            new = linspace(0, 1, length(B_ambig));
            colors_ambig = zeros(length(B_ambig), 3);
            for i=1:3
                colors_ambig(:,i) = min(max(interp1(old, base_ambig(:,i), new)', 0), 1);
            end
        elseif length(B_ambig)==3
            colors_ambig = base_ambig(:,:);
        elseif length(B_ambig)==2
            colors_ambig = base_ambig([1,3],:);
        else
            colors_ambig = base_ambig(1,:);
        end
    else
        colors_ambig = [];
    end
    colors = [colors_rel; colors_ambig];
    Colors = {};
    for c = 1:size(colors,1)
        Colors = [Colors colors(c,:)];
    end
    
    X = {};
    Y = {};
    for b = 1:length(B_rel)
        boundary = B_rel{b};
        x = new_fP(boundary(:,2));
        y = new_fA(boundary(:,1));
        x = [x x(1)];
        y = [y y(1)];
        X = [X x];
        Y = [Y y];
    end
    for b = 1:length(B_ambig)
        boundary = B_ambig{b};
        x = new_fP(boundary(:,2));
        y = new_fA(boundary(:,1));
        x = [x x(1)];
        y = [y y(1)];
        X = [X x];
        Y = [Y y];
    end

end