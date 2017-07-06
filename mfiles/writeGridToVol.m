function writeGridToVol(data, bbmin, bbmax, fileName)

fid = fopen(fileName, 'w');

a = 'VOL';
fwrite(fid, a(1), 'char');
fwrite(fid, a(2), 'char');
fwrite(fid, a(3), 'char');

fwrite(fid, char(3), 'char');
fwrite(fid, 1, 'int32'); % single-bit precision

fwrite(fid, size(data, 1), 'int32');
fwrite(fid, size(data, 2), 'int32');
fwrite(fid, size(data, 3), 'int32');

fwrite(fid, 1, 'int32'); % 1-channel

fwrite(fid, bbmin(1), 'single');
fwrite(fid, bbmin(2), 'single');
fwrite(fid, bbmin(3), 'single');

fwrite(fid, bbmax(1), 'single');
fwrite(fid, bbmax(2), 'single');
fwrite(fid, bbmax(3), 'single');

 
for k=1:size(data, 3)
    for j=1:size(data, 2)
        for i=1:size(data, 1)
            fwrite(fid, data(i, j, k), 'single');
        end
    end
end

fclose(fid);