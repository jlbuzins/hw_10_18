v1 = [1 2 3];
v2 = [5 6 7];
[dp, dd, ep, a, s] = vector_operate(v1,v2);

%case test 1
if (a == [6 8 10])
    fprintf('addition passed\n')
else
    fprintf('addition failed\n')
end

%case test 2
if (dp == [5 10 21])
    fprintf('dot product passed\n')
else
    fprintf('dot product failed\n')
end

function [dot_prod, dot_div, element_pow, add, sub] = vector_operate(x,y)
dot_prod = x.*y;
dot_div = x./y;
element_pow = x.^y;
add = x+y;
sub = x-y;
end