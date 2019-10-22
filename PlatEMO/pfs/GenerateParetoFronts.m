% Cian Steenkamp
% Adapted from https://github.com/iostream6/ParetoFrontGenerator
% The code in the link above is old and seems to not support the newest PlatEMO version.

clear;
clc;

dimensions = [3; 5; 8; 10; 15;]; % #objectives
paretoPoints = 1500; % not exactly so many points are generated AND it differs for each #obj

% ------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------
%    Generate DTLZ Pareto fronts
% ------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------

%% DTLZ1 pareto fronts
filenames = {"DTLZ1.3D.pf"; "DTLZ1.5D.pf"; "DTLZ1.8D.pf"; "DTLZ1.10D.pf"; "DTLZ1.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    PF = UniformPoint(paretoPoints, M) / 2;
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% DTLZ2 pareto fronts
filenames = {"DTLZ2.3D.pf"; "DTLZ2.5D.pf"; "DTLZ2.8D.pf"; "DTLZ2.10D.pf"; "DTLZ2.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    PF = UniformPoint(paretoPoints, M);
    PF = PF ./ repmat(sqrt(sum(PF .^ 2, 2)), 1, M);
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% DTLZ3 pareto fronts
filenames = {"DTLZ3.3D.pf"; "DTLZ3.5D.pf"; "DTLZ3.8D.pf"; "DTLZ3.10D.pf"; "DTLZ3.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    PF = UniformPoint(paretoPoints, M);
    PF = PF ./ repmat(sqrt(sum(PF .^ 2, 2)), 1, M);
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% DTLZ4 pareto fronts
filenames = {"DTLZ4.3D.pf"; "DTLZ4.5D.pf"; "DTLZ4.8D.pf"; "DTLZ4.10D.pf"; "DTLZ4.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    PF = UniformPoint(paretoPoints, M);
    PF = PF ./ repmat(sqrt(sum(PF .^ 2, 2)), 1, M);
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% DTLZ5 pareto fronts
filenames = {"DTLZ5.3D.pf"; "DTLZ5.5D.pf"; "DTLZ5.8D.pf"; "DTLZ5.10D.pf"; "DTLZ5.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    PF = [0:1 / (paretoPoints - 1):1; 1: - 1 / (paretoPoints - 1):0]';
    PF = PF ./ repmat(sqrt(sum(PF .^ 2, 2)), 1, size(PF, 2));
    PF = [PF(:, ones(1, M - 2)), PF];
    PF = PF ./ sqrt(2) .^ repmat([M - 2, M - 2: - 1:0], size(PF, 1), 1);
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% DTLZ6 pareto fronts
filenames = {"DTLZ6.3D.pf"; "DTLZ6.5D.pf"; "DTLZ6.8D.pf"; "DTLZ6.10D.pf"; "DTLZ6.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    PF = [0:1 / (paretoPoints - 1):1; 1: - 1 / (paretoPoints - 1):0]';
    PF = PF ./ repmat(sqrt(sum(PF .^ 2, 2)), 1, size(PF, 2));
    PF = [PF(:, ones(1, M - 2)), PF];
    PF = PF ./ sqrt(2) .^ repmat([M - 2, M - 2: - 1:0], size(PF, 1), 1);
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% DTLZ7 pareto fronts
filenames = {"DTLZ7.3D.pf"; "DTLZ7.5D.pf"; "DTLZ7.8D.pf"; "DTLZ7.10D.pf"; "DTLZ7.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    interval = [0, 0.251412, 0.631627, 0.859401];
    median = (interval(2) - interval(1)) / (interval(4) - interval(3) + interval(2) - interval(1));
    X = ReplicatePoint(paretoPoints, M - 1);
    X(X <= median) = X(X <= median) * (interval(2) - interval(1)) / median + interval(1);
    X(X > median) = (X(X > median) - median) * (interval(4) - interval(3)) / (1 - median) + interval(3);
    PF = [X, 2 * (M - sum(X / 2 .* (1 + sin(3 * pi .* X)), 2))];
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

% ------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------
%    Generate WFG Pareto fronts
% ------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------

%% WFG1 pareto fronts
filenames = {"WFG1.2D.pf"; "WFG1.3D.pf"; "WFG1.5D.pf"; "WFG1.10D.pf"; "WFG1.15D.pf"};

for z = 1:size(dimensions)
    M = dimensions(z);
    P = UniformPoint(paretoPoints, M);
    c = ones(size(P, 1), M);
    for i = 1 : size(P, 1)
        for j = 2 : M
            temp = P(i, j) / P(i, 1) * prod(1 - c(i, M - j + 2:M - 1));
            c(i, M - j + 1) = (temp ^ 2 - temp + sqrt(2 * temp)) / (temp ^ 2 + 1);
        end
    end
    x = acos(c) * 2 / pi;
    temp = (1 - sin(pi / 2 * x(:, 2))) .* P(:, M) ./ P(:, M - 1);
    a = 0 : 0.0001 : 1;
    E = abs(temp * (1 - cos(pi / 2 * a)) - 1 + repmat(a + cos(10 * pi * a + pi / 2) / 10 / pi, size(x, 1), 1));
    [~, rank] = sort(E, 2);
    for i = 1 : size(x, 1)
        x(i, 1) = a(min(rank(i, 1:10)));
    end
    P = convex(x);
    P(:, M) = mixed(x);
    PF = repmat(2:2:2 * M, size(P, 1), 1) .* P;
    dlmwrite(filenames{z}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% WFG2 pareto fronts
filenames = {"WFG2.2D.pf"; "WFG2.3D.pf"; "WFG2.5D.pf"; "WFG2.10D.pf"; "WFG2.15D.pf"};

for z = 1:size(dimensions)
    M = dimensions(z);
    P = UniformPoint(paretoPoints, M);
    c = ones(size(P, 1), M);
    for i = 1 : size(P, 1)
        for j = 2 : M
            temp = P(i, j) / P(i, 1) * prod(1 - c(i, M - j + 2:M - 1));
            c(i, M - j + 1) = (temp ^ 2 - temp + sqrt(2 * temp)) / (temp ^ 2 + 1);
        end
    end
    x = acos(c) * 2 / pi;
    temp = (1 - sin(pi / 2 * x(:, 2))) .* P(:, M) ./ P(:, M - 1);
    a = 0 : 0.0001 : 1;
    E = abs(temp * (1 - cos(pi / 2 * a)) - 1 + repmat(a .* cos(5 * pi * a) .^ 2, size(x, 1), 1));
    [~, rank] = sort(E, 2);
    for i = 1 : size(x, 1)
        x(i, 1) = a(min(rank(i, 1:10)));
    end
    P = convex(x);
    P(:, M) = disc(x);
    P = P(NDSort(P, 1) == 1, :);
    PF = repmat(2:2:2 * M, size(P, 1), 1) .* P;
    dlmwrite(filenames{z}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% WFG3 pareto fronts
filenames = {"WFG3.2D.pf"; "WFG3.3D.pf"; "WFG3.5D.pf"; "WFG3.10D.pf"; "WFG3.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    X = (0:1 / (paretoPoints - 1):1)';
    X = [X, zeros(paretoPoints, M - 2) + 0.5, zeros(paretoPoints, 1)];
    P = linear(X);
    PF = repmat(2:2:2 * M, size(P, 1), 1) .* P;
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% WFG4 pareto fronts
filenames = {"WFG4.2D.pf"; "WFG4.3D.pf"; "WFG4.5D.pf"; "WFG4.10D.pf"; "WFG4.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    P = UniformPoint(paretoPoints, M);
    P = P ./ repmat(sqrt(sum(P .^ 2, 2)), 1, M);
    PF = repmat(2:2:2 * M, size(P, 1), 1) .* P;
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% WFG5 pareto fronts
filenames = {"WFG5.2D.pf"; "WFG5.3D.pf"; "WFG5.5D.pf"; "WFG5.10D.pf"; "WFG5.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    P = UniformPoint(paretoPoints, M);
    P = P ./ repmat(sqrt(sum(P .^ 2, 2)), 1, M);
    PF = repmat(2:2:2 * M, size(P, 1), 1) .* P;
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% WFG6 pareto fronts
filenames = {"WFG6.2D.pf"; "WFG6.3D.pf"; "WFG6.5D.pf"; "WFG6.10D.pf"; "WFG6.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    P = UniformPoint(paretoPoints, M);
    P = P ./ repmat(sqrt(sum(P .^ 2, 2)), 1, M);
    PF = repmat(2:2:2 * M, size(P, 1), 1) .* P;
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% WFG7 pareto fronts
filenames = {"WFG7.2D.pf"; "WFG7.3D.pf"; "WFG7.5D.pf"; "WFG7.10D.pf"; "WFG7.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    P = UniformPoint(paretoPoints, M);
    P = P ./ repmat(sqrt(sum(P .^ 2, 2)), 1, M);
    PF = repmat(2:2:2 * M, size(P, 1), 1) .* P;
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% WFG8 pareto fronts
filenames = {"WFG8.2D.pf"; "WFG8.3D.pf"; "WFG8.5D.pf"; "WFG8.10D.pf"; "WFG8.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    P = UniformPoint(paretoPoints, M);
    P = P ./ repmat(sqrt(sum(P .^ 2, 2)), 1, M);
    PF = repmat(2:2:2 * M, size(P, 1), 1) .* P;
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

%% WFG9 pareto fronts
filenames = {"WFG9.2D.pf"; "WFG9.3D.pf"; "WFG9.5D.pf"; "WFG9.10D.pf"; "WFG9.15D.pf"};

for i = 1:size(dimensions)
    M = dimensions(i);
    P = UniformPoint(paretoPoints, M);
    P = P ./ repmat(sqrt(sum(P .^ 2, 2)), 1, M);
    PF = repmat(2:2:2 * M, size(P, 1), 1) .* P;
    dlmwrite(filenames{i}, PF, 'precision', '%.4f', 'delimiter', ' ');
end

function W = ReplicatePoint(SampleNum,M)
    if M > 1
        SampleNum = (ceil(SampleNum^(1/M)))^M;
        Gap       = 0:1/(SampleNum^(1/M)-1):1;
        eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
        eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
    else
        W = (0:1/(SampleNum-1):1)';
    end
end

function Output = convex(x)
    Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = mixed(x)
    Output = 1-x(:,1)-cos(10*pi*x(:,1)+pi/2)/10/pi;
end

function Output = disc(x)
    Output = 1-x(:,1).*(cos(5*pi*x(:,1))).^2;
end

function Output = linear(x)
    Output = fliplr(cumprod([ones(size(x,1),1),x(:,1:end-1)],2)).*[ones(size(x,1),1),1-x(:,end-1:-1:1)];
end

