function P = my_reset_orientation(P,MM)

if ~isempty(P(1).private.extras) && isstruct(P(1).private.extras) && isfield(P(1).private.extras,'mat')
    for i=1:size(P(1).private.dat,4)
        mat(:,:,i) = MM;
    end
end
for k=1:numel(P)
    P(k).mat = MM;
    P(k).private.mat = MM;
    if ~isempty(P(k).private.extras) && isstruct(P(k).private.extras) && isfield(P(k).private.extras,'mat')
        P(k).private.extras.mat = mat;
    end
end