function M = vector2skewsym(x)
% Find w^ from w

    if size(x,2)>1 
        disp('error: use column vector');
        return
    elseif size(x,1)~=3
        disp('error: vector has length not equal to 3');
        return
    else
        M = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
    end

end
