function vol = load3DTif_uint16(path)
%From Dan Goodwin's ExSeqProcessing repository: https://github.com/dgoodwin208/ExSeqProcessing/blob/master/utils/load3DImage_uint16.m
    vol = double(read_file(path));
end
