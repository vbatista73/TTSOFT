function [y]=cellSplit(x,ds)
% This subroutine breaks up a TT-tensor into subtensors corresponding to
% each physical dimension and stores them in a cell
D = numel(ds);
y = {};
for i=1:D
    y{i} = chunk(x,1,ds(i));
    if(i<D); x=chunk(x,ds(i)+1,x.d); end;
end

end