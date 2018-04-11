%% Network structure (HOMEWORK 1)
%% Graph generation

clear
clc
load routes.mat

%create EdgeTable without duplicates
EdgeTable = unique(table([(routes.sourceairport) (routes.destinationairport)],... 
    'VariableNames', {'EndNodes'}));
%add weights to EdgeTable using concatenation
numedges = size(EdgeTable,1);
EdgeTable = [EdgeTable table(ones(numedges,1),'VariableNames',{'Weight'})];
EdgeTable.EndNodes = cellstr(EdgeTable.EndNodes);

%create directed graph
G = digraph(EdgeTable);
clear EdgeTable numedges routes

%% Plot graph

figure;
plot1 = plot(G,'Layout', 'force');
plot1.NodeColor = 'k';
plot1.EdgeColor = 'r';
plot1.LineStyle = '--';
title('Fig.1: Airport network');
%plot.NodeLabel = G.Nodes.Name;     %uncomment to display labels

%% Network properties and parameters

[H, GC, Hubs] = networkProperties(G);

%% Network communities (HOMEWORK 2)
%% Link analysis

[rank, iterations] = pagerank(G);
ScoresTable = table(G.Nodes.Name, rank, centrality(G,'pagerank'),...
    centrality(G,'hubs'), centrality(G,'authorities'),...
    'VariableNames', {'Node', 'Rank_Man', 'Rank', 'Hub', 'Authority'});
clear rank

%plot
figure
hold on
scatter(1:numnodes(G), ScoresTable.Rank);
scatter(1:numnodes(G), ScoresTable.Rank_Man, 'xk');
title('Fig.1: ranks');
ylabel('rank');
xlabel('node index');
legend('PageRank', 'Rank');
xlim([0,numnodes(G)+1]);
hold off

figure
hold on
scatter(1:numnodes(G), ScoresTable.Hub);
scatter(1:numnodes(G), ScoresTable.Authority, 'xk');
title('Fig.2: hubs and authorities');
ylabel('scores');
xlabel('node index');
legend('Hub score', 'Authority score');
xlim([0,numnodes(G)+1]);
hold off


%top 16 nodes (by score) (following values are set to get 16 elements)
top_hub = G.Nodes.Name(ScoresTable.Hub > 0.0053);
top_authority = G.Nodes.Name(ScoresTable.Authority > 0.0052);
top_rank_man = G.Nodes.Name(ScoresTable.Rank_Man > 0.003);
top_rank_aut = G.Nodes.Name(ScoresTable.Rank > 0.003);
TopRank = table(top_rank_man, top_rank_aut, top_hub, top_authority,...
    'VariableNames', {'Rank_Manual','Rank_MATLAB','Hub_MATLAB','Authority_MATLAB'});
clear top_rank_man top_rank_aut top_hub top_authority
fprintf('The algorithm pagerank converges in %d iterations\n', iterations)
TopRank;

%% Link prediction

S = similarity_PA(G);

%% Community detection

%compute only on giant component because other components are already separated
Clusters = spectralclustering(GC);
cluster1 = Clusters.Name(Clusters.Cluster == 1);
cluster2 = Clusters.Name(Clusters.Cluster == 2);
cluster3 = Clusters.Name(Clusters.Cluster == 3);

%plot GC
figure
plot2 = plot(GC,'Layout','force');
plot2.LineStyle = '--';
highlight(plot2, findnode(GC, cluster1), 'NodeColor', 'k')
highlight(plot2, findnode(GC, cluster2), 'NodeColor', 'g')
highlight(plot2, findnode(GC, cluster3), 'NodeColor', 'r')
title('Fig.5: airport network split in communities');

