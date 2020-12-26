function A = gramschmidt(M,V)
% Construct an orthogonal basis in the space defined by M
A = zeros(size(M,1),size(M,2));
A(:,1) = V;
A(:,1) = A(:,1)/norm(A(:,1),2);

for i = 2:size(A,2)
    tmp = zeros(size(M,1),1);
    for j = 1:i-1
        tmp = tmp + A(:,j)'*M(:,i)/(A(:,j)'*A(:,j))*A(:,j);
    end
        A(:,i) = M(:,i) - tmp;
        A(:,i) = A(:,i)/norm(A(:,i),2);
end
end



