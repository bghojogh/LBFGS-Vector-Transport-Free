function arr = convert_struct_to_array(struct_, field_)
    i = 1;
    arr = [];
    while 1
        try
            arr(i) = getfield(struct_(i),field_);
            i = i + 1;
        catch
            break
        end
    end
end