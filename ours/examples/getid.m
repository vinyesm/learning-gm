function x=getid(cc)
    c=char(cc);
    cs=strsplit(c,',');
    x=str2num(cs{1});
%     keyboard;
end