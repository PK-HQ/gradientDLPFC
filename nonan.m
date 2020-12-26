function nonanmat=nonan(mat)
%RM NAN
nonanmat = mat(~isnan(mat));
end