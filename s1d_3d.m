function y=s1d_3d(x,remove0)
if nargin<2
remove0=1;
end
    if remove0==1
    x(find(x==0))=[];
    y=x;
    y=x(:);
    else
     y=x(:);
    end
    