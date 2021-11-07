function writeELEM(ELEM, filename)

if nargin<2
    filename = 'ELEM.dat';
end

% write ELEM array into a file
nel = size(ELEM, 2); % number of elements
nn = size(ELEM, 1); % node per element

fid = fopen(filename, 'w');

fprintf(fid, '%d %d\n',  [nel, nn]);

% create format for connectivity
fmt = '';
for i = 1:nn
    fmt = [fmt,'%d '];
end

fmt = strip(fmt);
fmt = [fmt,'\n'];    
fprintf(fid, fmt, ELEM);
fclose(fid);
end

