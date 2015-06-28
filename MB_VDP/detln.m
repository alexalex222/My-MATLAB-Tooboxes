function [Y, error_flag] = detln( X )
% Y = logdet( X )
% return a log determinant of X
[d err] = chol(X);
if err
  error('error in Choleski disposition for logdet')
  return
end
d = diag(d);
Y = sum( log(d) ) *2;
error_flag=0;

%%%%%% Local Variables: *****
%%%%%% mode: matlab *****
%%%%%% End: *****
