function [n,nL,nR,nT,nB] = Final_Index(N,i,j)
    n = i + (j - 1) * N;

    if i > 1
        nL = i - 1 + (j - 1) * N;
    else
        nL = i + N - 2 + (j - 1) * N;
    end

    if i < N
        nR = i + 1 + (j - 1) * N;
    else
        nR = i - N + 2 + (j - 1) * N;
    end

    if j < N
        nT = i + j * N;
    else
        nT = i + (j - N + 2) * N;
    end

    if j > 1
        nB = i + (j - 2) * N;
    else
        nB = i + (j + N - 2) * N;
    end
end