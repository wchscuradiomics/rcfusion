function p = fisher_test_R_C(x, num)
%  Input variable: x
    %   x is the contingency table. You can use function 'crosstab'  to generate the contingency table.
    %   num : the number of Monte Carlo simulations to be performed, default is 100000.
%  Output variable: p
    %    p：p-values
% Source of subfunction 'rcont2':  https://people.math.sc.edu/Burkardt/m_src/asa144/asa144.html
% for example
% x = [    11    11
%      0     6
%     11    17];
% p = fisher_test_R_C(x)  
% p =  0.0977;
% Because of the randomness in the algorithm, the result may be different each time.
% You can increase the value of  'num'  to make the result more robust.
    if nargin == 1
        num = 100000;
    end
    [r,c] = size(x);
    m = sum(x,2);
    n =  sum(x,1);
    nrow= r;
    ncol= c;
    nrowt = m;
    ncolt = n;
    xxx = prod(prod(factorial(x)));
    rcont2( nrow, ncol, nrowt, ncolt,false,randi(10000000));   % Initialization, cannot be deleted
    js = 0;
    for ii = 1:num
        [ ~, ~, y] = rcont2( nrow, ncol, nrowt, ncolt,true,randi(10000000));
        yyy = prod(prod(factorial(y)));
        if yyy >= xxx
            js = js+1;
        end
    end
    p = js/num;
end