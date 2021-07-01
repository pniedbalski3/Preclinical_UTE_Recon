function line_out = check_method_line_read(line_in)

%New Bruker convention is to use @*(0) to indicate a bunch of zeros...
%Matlab doesn't handle this, so need to check

if contains(line_in,'@')
    ind = strfind(line_in,'@');
    ind_end1 = strfind(line_in,'*');
    ind_end2 = strfind(line_in,')');
    n_repeat = str2num(line_in((ind+1):(ind_end1-1)));
    val = str2num(line_in((ind_end1+1):ind_end2));
    replace_num = ones(1,n_repeat)*val;
    replace_str = num2str(replace_num);
    old_str = line_in(ind:ind_end2);
    line_out = strrep(line_in,old_str,replace_str);
else
    line_out = line_in;
end

