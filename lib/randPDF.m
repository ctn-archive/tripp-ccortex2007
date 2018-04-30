% function result = randPDF(M, N, PDF, step) generates an M-by-N matrix with random
% entries drawn from the discrete probability density function 'PDF' with
% step size 'step'. 

function result = randPDF(M, N, PDF, step) 
    result = [];
    c = cumsum(PDF);
    c = c / c(end);
    r = rand(M, N);
    for i = 1:M
        for j = 1:N   
            block = max([0, find(c < r(i,j))]) + 1;
            result(i,j) = step * block;
        end
    end
    