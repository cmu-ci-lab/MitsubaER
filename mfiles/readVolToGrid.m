function [data, bbmin, bbmax] = readVolToGrid(fileName)
data = 0;
bbmin = zeros(3, 1);
bbmax = zeros(3, 1);
fid = fopen(fileName, 'r');

header = char(3);
header(1) = fread(fid, 1, 'char');
header(2) = fread(fid, 1, 'char');
header(3) = fread(fid, 1, 'char');

if(~strcmp(header, 'VOL'))
    disp('failed to read header');
    return;
end

version = fread(fid, 1, 'char');
if(version ~= 3)
    disp('wrong version');
    return;
end
type = fread(fid, 1, 'int32'); % single-bit precision
if(type ~= 1)
    disp('can only handle single precision');
    return;
end

xres = fread(fid, 1, 'int32');
yres = fread(fid, 1, 'int32');
zres = fread(fid, 1, 'int32');

data = zeros(xres, yres, zres);

channels = fread(fid, 1, 'int32'); % 1-channel
if(channels ~= 1)
    disp('can only single channel');
    return;
end

bbmin(1) = fread(fid, 1, 'single');
bbmin(2) = fread(fid, 1, 'single');
bbmin(3) = fread(fid, 1, 'single');

bbmax(1) = fread(fid, 1, 'single');
bbmax(2) = fread(fid, 1, 'single');
bbmax(3) = fread(fid, 1, 'single');

for k=1:zres
    for j=1:yres
        for i=1:xres
            data(i, j, k) = fread(fid, 1, 'single');
        end
    end
end

fclose(fid);