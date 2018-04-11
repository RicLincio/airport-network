function Clusters = spectralclustering(J)

A = adjacency(J); A = A';
D = diag(sum(A,2));
T = inv(D.^(0.5));
L = sparse(eye(numnodes(J)) - T * A * T);
[X,E] = eigs(L,30,'smallestabs');    %heuristically k = 3
e = diag(E);
figure
scatter(1:30,e);
title('Fig.3: the smallest eigenvalues of L');
xlabel('index');
ylabel('value');
V = T * X;

%modify for greater K
k = 2;      %excluding the zero eigenvalue
V = V(:,2:k+1);
for i = 1:k
    V(:,i) = V(:,i) / norm(V(:,i));
end
figure
scatter(V(:,1), V(:,2));
title('Fig.4: Bidimensional representation of nodes');
xlim([-0.02 0.2])
xlabel('v_{N-1}');
ylabel('v_{N-2}');


[C, mu] = kmeans(V,3,k);	%from the scatter we can distinguish 3 main clusters
Clusters = table(J.Nodes.Name, C, 'VariableNames', {'Name','Cluster'});

end