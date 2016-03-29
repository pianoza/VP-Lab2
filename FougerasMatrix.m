% Get transformation matrix using Fougeras method.
function [ fougerasMatrix, intrinsics, extrinsics ] = FougerasMatrix( points3D, points2D )
    pN = length(points3D);
    Q = zeros(2*pN, 11);
    B = zeros(2*pN, 1);
    for i = 1:pN
        Q(2*i-1, :) = [points3D(i,:), -points2D(1,i)*points3D(i,:), 0., 0., 0., 1., 0.];
        Q(2*i, :) = [0, 0, 0, -points2D(2,i)*points3D(i,:), points3D(i,:), 0., 1.];
        B(2*i-1) = points2D(1,i);
        B(2*i) = points2D(2,i);
    end;
    X = Q\B;
    T1 = X(1:3)';
    T2 = X(4:6)';
    T3 = X(7:9)';
    C1 = X(10);
    C2 = X(11);
    %compute intrinsic parameters:
    normT2 = norm(T2,2)^2;
    U0 = (T1*T2')/normT2;
    V0 = (T2*T3')/normT2;
    Au = norm(cross(T1',T2'))/normT2;
    Av = norm(cross(T2',T3'))/normT2;
    intrinsics = [Au, 0, U0, 0; 0, Av, V0, 0; 0, 0, 1, 0];
    % intMat and intMat2 are the same;
    %compute extrinsic parameters:
    r1 = (norm(T2)/norm(cross(T1',T2')))*(T1-(T1*T2'/normT2)*T2);
    r2 = (norm(T2)/norm(cross(T2',T3')))*(T3-(T2*T3'/normT2)*T2);
    r3 = T2/norm(T2);
    Tx = (norm(T2)/norm(cross(T1',T2')))*(C1-(T1*T2'/normT2));
    Ty = (norm(T2)/norm(cross(T2',T3')))*(C2-(T2*T3'/normT2));
    Tz = 1/norm(T2);
    extrinsics = [
        r1, Tx;
        r2, Ty;
        r3, Tz;
        0,0,0,1
        ];
    fougerasMatrix = intrinsics*extrinsics;
    fougerasMatrix = fougerasMatrix./fougerasMatrix(3, 4);
end

