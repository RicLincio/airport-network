% Network properties and parameters

function [H, GC, Hubs] = networkProperties(G)

%sparsity
if numedges(G) < numnodes(G) *(numnodes(G)-1) / 100
    display('The network is sparse');
else
    display('The network is not sparse');
end

%degree of nodes
G.Nodes.Degree = indegree(G) + outdegree(G);

%average degree
fprintf("The average degree is %3.2f \n", mean(G.Nodes.Degree))

%degree variance
fprintf("The degree variance is %4.2f \n", var(G.Nodes.Degree))

%check presence of isolated nodes
isolatedNodes = 0;
for j = 1:numnodes(G)
    if G.Nodes.Degree(j) == 0
        isolatedNodes = isolatedNodes + 1;
    end
end
fprintf("There are %i isolated nodes\n", isolatedNodes)
clear j isolatedNodes

%hub
[k_hub, hub_index] = max(G.Nodes.Degree);
fprintf("The biggest hub is %s with degree %i \n", string(G.Nodes.Name(hub_index)), k_hub)
hd = 300;       %minimum degree value to be considered a hub
Hubs = table(G.Nodes.Name(G.Nodes.Degree > hd), G.Nodes.Degree(G.Nodes.Degree > hd),...
    'VariableNames', {'Name', 'Degree'});
Hubs = sortrows(Hubs, 'Degree', 'descend');
clear hd k_hub hub_index

%degree distribution
p = zeros(Hubs.Degree(1), 1);  %max degree is 2*(N-1) for directed graphs
for j = 1:numnodes(G)
    p(G.Nodes.Degree(j)) = p(G.Nodes.Degree(j)) + 1;
end
p = p / numnodes(G);
clear j k

%plot degree distribution
figure
%scatter(k_log, p_log, 'x')
scatter(1:1:length(p), p, 'x');
set(gca, 'xscale', 'log', 'yscale', 'log')
xlim([1,1000])
grid on
title('Fig.4: Degree distribution')
xlabel('degree (k)')
ylabel('frequency p(k)')

%logarithmic binning
W = 2;
p_log = zeros(9,1); %set 9 as length because (W^9 = 2^9 > k_hub)
for j = 1:length(p_log)
    if W^j < length(p)
        p_log(j) = sum(p(W^(j-1):W^(j)-1));
    else
        p_log(j) = sum(p(W^(j-1):length(p)));
    end
end
y = zeros(1,length(p)); %recentering
for j = 1:length(p_log)
    index = floor(W^(j-1) + ((W^j - W^(j-1)) / 2));
    y(index) = p_log(j);
end
figure 
scatter(1:1:length(p), y/numnodes(G), 'x')
set(gca, 'xscale', 'log', 'yscale', 'log')
xlim([1,1000])
grid on
title('Fig.5: Degree distribution with logarithmic binning')
xlabel('degree (k)')
ylabel('frequency p(k)')
%hold on
clear W y j index p p_log

%gamma coefficient
gamma_approx = 1 + numnodes(G) / (sum(log(G.Nodes.Degree / min(G.Nodes.Degree))));
fprintf('The approximated value of gamma is %1.4f\n', gamma_approx)
%gamma_estim
clear gamma_approx

%neighbors
%nearest(G,'VCE',1);

%connectivity
%span the giant component and create subgraph
GC = bfsearch(G, Hubs.Name(1));
fprintf("The giant component has %i nodes (%1.3f%%)\n", length(GC(:,1)), length(GC(:,1))/numnodes(G))
GC = subgraph(G, GC);

%reorder nodes with GC first
order = bfsearch(G, Hubs.Name(1), 'Restart', true);
H = reordernodes(G,order);
imshow(abs(full(adjacency(H)) - 1));
title('Fig.3: Adjacency matrix')
ylabel('y: departure airport')
xlabel('x: arrival airport')
clear order

%find number of connected components (with strongly connected nodes)
bins = conncomp(G);
components = 0;
for j = 1:length(bins)
    if bins(j) > components
        components = bins(j);
    end
end
fprintf("The number of connected components (with strongly connected nodes) is %i \n", components)
clear bins components j

%small world properties
d_G = distances(G);
d_G = d_G(isfinite(d_G));       %remove elements with value infinite
%giant component
d_GC = distances(GC);
d_GC = d_GC(isfinite(d_GC));	%remove elements with value infinite
fprintf('The expected average distance of G is %1.3f\n', log(numnodes(G)) / log(mean(G.Nodes.Degree)))
fprintf('The average distance of G is %1.3f\n', mean(d_G))
fprintf('The average distance of GC is %1.3f\n', mean(d_GC))
fprintf('The diameter of G is %i\n', max(d_G))
fprintf('The diameter of GC is %i \n', max(d_GC))
clear d_G d_GC

%clustering coefficient
for j = 1:numnodes(G)     %for each node
    neighbor = unique([predecessors(G, j); successors(G, j)]);
    NB = subgraph(G, neighbor);
    G.Nodes.ClusteringCoeff(j) = numedges(NB) / (numnodes(NB) * (numnodes(NB) - 1));  %compute clustering coefficient (directed case)
end
C_finite = G.Nodes.ClusteringCoeff(isfinite(G.Nodes.ClusteringCoeff)); %compute only on finite elements
fprintf('The average clustering coefficient is %1.3f \n', mean(C_finite))
clear neighbor j C_finite NB

end