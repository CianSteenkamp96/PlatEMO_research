function [FrontNo,MaxFNo] = NDSort_CDAS(PopObj,nSort)

% My understanding: the archive = the first non-dominated front (NDSort does non-dominated sort on a set of individuals - resulting in different non-dominated front layers). That is, the max size of the
% arhive will be the same as the population size. See NDSort.m for more
% info.

% Copyright 2016-2017 Ye Tian

    [N,M] = size(PopObj);

    %% Parameter setting
    S = 0.39;

    %% Transformation of objective values
    Norm   = sqrt(sum(PopObj.^2,2));
    Angle  = acos(1-pdist2(PopObj,eye(M),'cosine'));
    PopObj = repmat(Norm,1,M).*sin(Angle+S*pi)./sin(S*pi);
    [PopObj,rank] = sortrows(PopObj);
    
    %% Non-dominated sorting
    FrontNo = inf(1,N);
    MaxFNo  = 0;
    while sum(FrontNo<inf) < min(nSort,N)
        MaxFNo = MaxFNo + 1;
        for i = 1 : N
            if FrontNo(i) == inf
                Dominated = false;
                for j = i-1 : -1 : 1
                    if FrontNo(j) == MaxFNo
                        m = 2;
                        while m <= M && PopObj(i,m) >= PopObj(j,m)
                            m = m + 1;
                        end
                        Dominated = m > M;
                        if Dominated || M == 2
                            break;
                        end
                    end
                end
                if ~Dominated
                    FrontNo(i) = MaxFNo;
                end
            end
        end
    end
    FrontNo(rank) = FrontNo;
end
