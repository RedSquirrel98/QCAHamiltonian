function tensor= tensor_product(index,N,sigma)


    if index==1
        tensor=sigma;
    elseif index >=2
        tensor=kron(eye(2^(index-1)),sigma);
    end


    rem=N-index;
    tensor=kron(tensor,eye(2^rem));


end

