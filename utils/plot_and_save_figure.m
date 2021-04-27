function plot_and_save_figure(arr, name_, y_label, path_)

if ~exist(path_, 'dir')
   mkdir(path_)
end

f = figure('visible', 'off');
plot(arr, "-*")
ylabel(y_label); xlabel("iteration");
print(path_+name_+'.png', '-dpng', '-r100');
close(f)

end