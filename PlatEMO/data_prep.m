

a = importdata('NSGAIII_DTLZ7_M3_D1000_3.mat');

numObj = length(a.result{2}(1).obj);
popSize = length(a.result{2});

x = zeros(popSize, numObj);
for i = 1:popSize
    for j = 1:numObj
        x(i, j) = a.result{2}(i).obj(j);
    end    
end 

dlmwrite('NSGAIII_DTLZ7_M3_D1000_3_obj.out', x, ' ');