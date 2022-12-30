function t = twistt(v)
% Find twist for translation case
    
    if size(v,2)>1
        disp('error: use row vector');
        return
    elseif size(v,1)~=3
        disp('error: vector has length not equal to 3');
        return
    else
        % Find the twist
        t=[v; 0; 0; 0];
    end
    
end
