function outMat = householder(inMat)
%  Performs a Householder transformation on the input matrix H
%     - Saves upper-triangularized H as H for output
%     - Takes advantage of Matlab's matrix multiplication routines

% Initialization
outMat = inMat;
[m, n] = size(outMat);

% Iterate through each column up to the min(m, n)
for k = 1:min(m, n) 
    % Extract the vector to be transformed
    z = outMat(k:m, k);
    sigma = norm(z) * sign(z(1));
    
    if sigma ~= 0
        u = z;
        u(1) = u(1) + sigma;
        u = u / norm(u);
        
        % Apply Householder transformation
        outMat(k:m, k:n) = outMat(k:m, k:n) - 2 * (u * (u' * outMat(k:m, k:n)));
        
        % Ensure numerical stability by enforcing zeros explicitly
        outMat(k+1:m, k) = 0;
    end
end
end
