% Get Hall transformation matrix
function hallMatrix = HallMatrix(points3D, points2D)
    pN = length(points3D);
    Q = zeros(2*pN, 11);
    B = zeros(2*pN, 1);
    for i = 1 : pN
        Q(2*i-1,:)=[points3D(i, 1), points3D(i, 2), points3D(i, 3), 1, 0, 0, 0, 0, -points2D(1,i)*points3D(i,1), -points2D(1,i)*points3D(i,2), -points2D(1,i)*points3D(i,3)];
        Q(2*i,:)=[0, 0, 0, 0, points3D(i, 1), points3D(i, 2), points3D(i, 3), 1, -points2D(2,i)*points3D(i,1), -points2D(2,i)*points3D(i,2), -points2D(2,i)*points3D(i,3)];
        B(2*i-1) = points2D(1, i);
        B(2*i) = points2D(2, i);
    end;
    A = Q\B;
    A(12) = 1.0;
    A = reshape(A, [4, 3]);
    hallMatrix = A';
end