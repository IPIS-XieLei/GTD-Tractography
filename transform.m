function Tracts = transform(streamlines, affine)

Tracts = cell(1,1);

% Transpose
if size(streamlines,1)==1
    streamlines = streamlines';
end

% Start convert
for i = 1:size(streamlines,1)
    Tract = streamlines(i,:);
    t = Tract{1}; t(:,4) = 1;
    Tract = t*affine';
    Tracts{i,1} = Tract(:,1:3);
end

end