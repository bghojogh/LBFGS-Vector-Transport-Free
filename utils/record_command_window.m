function record_command_window(path_save, start)
    if ~exist(path_save, 'dir')
       mkdir(path_save)
    end
    if start == "on"
        fid = fopen(path_save+"command_window.txt", 'at');
        fprintf(fid, '======================================================');
        fclose(fid);
        diary(path_save+"command_window.txt");  %--> save command window --> https://stackoverflow.com/questions/5833356/how-to-save-the-contents-of-matlabs-command-window-to-a-file
        diary('on');  %--> recording the command window
    elseif start == "off"
        diary('off');
    end
end