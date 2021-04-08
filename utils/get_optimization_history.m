function [cost_list, grad_norm_list, stepsize_list, time_list, time_iterations] = get_optimization_history(info_)

i = 1;
cost_list = [];
grad_norm_list = [];
stepsize_list = [];
time_list = [];
while 1
    try
        cost_list(i) = getfield(info_(i),"cost");
        grad_norm_list(i) = getfield(info_(i),"gradnorm");
        stepsize_list(i) = getfield(info_(i),"stepsize");
        time_list(i) = getfield(info_(i),"time");
        i = i + 1;
    catch
        break
    end
end
time_iterations = [time_list(2) - time_list(1)];
for itr = 2:length(time_list)
    time_iterations(itr) = time_list(itr) - time_list(itr-1);
end

end