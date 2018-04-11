%% log binning

k_log = log(1:1:length(p));
k_min = 3;
k_log_2 = k_log(k_min:end);
X = ones(length(k_log_2), 2);
X(:,2) = k_log_2;
%Y = 


%% gamma

k_min = min(k);
gamma = 1 + N / (sum(log(k / k_min)));

% pow_val = (-gamma * (1:1:length(p)));
% figure
% plot((1:1:length(p)), pow_val)
% scatter(1:1:length(p), pow_val, 'x')
% set(gca, 'xscale', 'log', 'yscale', 'log')
% %xlim([1,1000])
% grid on
% title('Degree distribution with logarithmic binning')
% xlabel('degree (k)')
% ylabel('frequency p(k)')