i = 1;
cost_list = [];
while 1
    try
        cost_list(i) = getfield(results(i),"cost");
        i = i + 1;
    catch
        break
    end
end
plot(cost_list);