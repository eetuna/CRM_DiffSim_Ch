function R = so3rot(w,t)
% Map so(3) to SO(3)
% From eq 2.14 in [1].

    % If norm-2 of w is less than eps then we will return an identity matrix
    eps = 1e-16;

    % Check dimension of w
    if size(w,2)>1
        disp('error: use row vector');
        return
    elseif size(w,1)~=3
        disp('error: vector has length not equal to 3');
        return
    else
        % If w smaller than eps
        if norm(w * t, 2) < eps
            R = eye(3);
        % Otherwise calculate rotation matrix
        else
            % Find w^
            w1 = vector2skewsym(w);
            % Comput R using eq 2.14
            R = eye(3) + w1/sqrt(w'*w)*sin(sqrt(w'*w)*t)+w1*w1/sqrt(w'*w)^2*(1-cos(sqrt(w'*w)*t));
        end
    end
    
end
