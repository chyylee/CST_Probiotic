function dy = gLV_jacobian(t,y,params)
    % Defines number of equations based on y
    neq = length(y);
    dy = zeros(neq,1);

    % Reformat parameters to easily be able to loop through
    p_grow = params(1:neq);
    p_int = reshape(params(neq+1:end),[neq neq]);

    dy = zeros(neq,neq);
    temp = zeros(neq,1); 
    for i  = 1:neq
        for j =  1:neq
            if i ~= j
                dy(i,j) = p_int(i,j)*y(i);
                temp(i) = temp(i) + p_int(i,j)*y(j); % save jacobian pattern for diag
            else
                dy(i,j) = p_grow(i) + 2*p_int(i,j)*y(i); % diag pattern
            end    
        end
    end
    dy = dy + diag(temp); 
end