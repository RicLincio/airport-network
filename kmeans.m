function [out, mu] = kmeans(V, k, dim)

%initialize centroids (3 points, one for each cluster) (no random init)
mu = [V(1,:);
      V(3366,:);
      V(3378,:)];

%iteration step
SAD = 1e10;
convergence = false;
while ~(convergence)
    C = zeros(length(V),2);
    for i = 1:length(V)             %for each point
        d = zeros(k,1);             %compute distance from each centroid
        for j = 1:k
            d(j) = norm(V(i,:) - mu(j,:));
        end
        [C(i,2),C(i,1)] = min(d);   %and keep the minimum
    end
    d_tot = sum(C(:,2));
    
    if (SAD - d_tot) / d_tot < 1e-3 %converges
        convergence = true;
        
    else                            %update centroids and SAD
        mu = zeros(k,dim);
        count = zeros(k,1);
        for i = 1:length(C)
            mu(C(i,1),:) = mu(C(i,1),:) + V(i,:);
            count(C(i,1)) = count(C(i,1)) + 1;
        end
        for i = 1:k
            mu(i,:) = mu(i,:) / count(i);
        end
        SAD = d_tot;
    end
end
out = C(:,1);
end