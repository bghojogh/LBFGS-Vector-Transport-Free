function sim1_gen_data(DIM, K, N, e, RESFOLDER)

if nargin < 4
    e = 10; % eccentricity
end

if nargin < 5
    RESFOLDER = 'result';
end

if nargin==0
    DIM = 2; % data dimensions
end


% degrees of separation between the cluster means
SEP(1).name = 'no';
SEP(1).factor = 0;
SEP(2).name = 'mid';
SEP(2).factor = 1;
SEP(3).name = 'high';
SEP(3).factor = 5;
SEP(4).name = 'low';
SEP(4).factor = 0.2;

for i = 1:numel(SEP)
    
    [X,ignore,ignore,gen_ll,theta_gen] = mixgen(N,N,K,DIM,SEP(i).factor,e);
    % save
    data = X(1:N , :).';
    sep = SEP(i).name;
    filename = sprintf('%s/sim1_data_dim(%d)_K(%d)_N(%d)_E(%d)_sep(%s)', RESFOLDER, DIM, K, N, e, sep);
    save(filename, 'data', 'gen_ll', 'sep', 'DIM', 'K', 'theta_gen')

    
    figure('Name',filename)
    if DIM >= 3
        scatter3(data(1,:), data(2,:), data(3,:))
    else
        scatter(data(1,:), data(2,:))
    end
    drawnow
end
