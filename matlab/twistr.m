function t = twistr(w,q)
% Find twist for rotation case
    
    if size(w,2)>1 || size(q,2)>1
        disp('error: use row vector');
        return
    elseif size(w,1)~=3 || size(q,1)~=3
        disp('error: vector has length not equal to 3');
        return
    else
        % Find the twist
        t=[-cross(w,q);w];
    end
    
end
