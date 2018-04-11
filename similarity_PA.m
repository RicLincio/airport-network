function S = linkpredict_PA(J)

k = outdegree(J);
A = adjacency(J);
S = zeros(numnodes(J), numnodes(J));
for i = 1:numnodes(J)
    for j = (i+1):numnodes(J)
        if A(i,j) == 0
            S(i,j) = k(i) * k(j);
        end
    end
end
S = S / max(max(S));
S = S + S'; %the result is symmetric

end