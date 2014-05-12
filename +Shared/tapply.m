function r = tapply( x )
%
% function r = tapply( x )
%
% Mimics R's tapply().
%
% x has 2 columns: a list of factors in column1, and values in column2.
% 
% DHO.
%
x = sortrows([x(:,1) x(:,2)],1);
u = unique(x(:,1));
r = transpose(1:size(unique(x(:,1)),1));
k = 1;
for i=1:size(u)
    tmp=u(i);
    f = find(x(:,1)==tmp);
    r(k) = mean(x(f,2));
    k = k+1;
end

r = [u r];