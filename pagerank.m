function [rank, iterations] = pagerank(J)

%initialization
t = 0;          %iteration number
epsilon = 1e-15;   %convergence condition
p_prev = zeros(numnodes(J),1);
p_curr = ones(numnodes(J),1) / numnodes(J);

%page rank with restarts
A = full(adjacency(J)); A = A';
q_0 = ones(numnodes(J),1);  %equal probability to jump from dead ends to any other node
e = outdegree(J); e = double(e == 0);
A = A + q_0 * e';
d = sum(A);
for j = 1:numnodes(J)
	M(:,j) = A(:,j) / d(j);
end
q_1 = ones(numnodes(J),1) / numnodes(J);
M = 0.85 * M + 0.15 * q_1 * ones(1,numnodes(J));

%iteration
while norm(p_curr - p_prev) > epsilon
    t = t + 1;
    p_prev = p_curr;
    p_curr = M * p_prev;
end

%output
iterations = t;
rank = p_curr;

end