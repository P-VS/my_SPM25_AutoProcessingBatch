function [lcm] = compute_lcm(x,y)
%COMPUTE_LMC(X,Y) Computes the least common denominator 
    
    x = int32(x);
    y = int32(y);
    %choose the greter number
    if x > y
        greater = x;
    else
        greater = y;
    end

    while (true)
        if (mod(greater,x) == 0 && mod(greater,y) == 0)

            lcm = greater;
            break;
        end
        greater = greater + 1;
    end

end

