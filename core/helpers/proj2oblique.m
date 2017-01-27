function [ Anorm ] = proj2oblique( A )
%PROJ2OBLIQUE   Normalize each slice to lie on the sphere.

    Anorm = NaN(size(A));
    for i = 1:size(A,3)
        Anorm(:,:,i) = A(:,:,i)/norm(A(:,:,i),'fro');
    end

end

