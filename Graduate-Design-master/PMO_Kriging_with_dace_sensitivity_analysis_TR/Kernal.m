function Kernal( M,F,n)
n = 4;
if M == []
    return F;
else
	if length(M) == 1
        a = M(1);
        if HDMR(a) ~= []
            return HDMR(a);
        else
            k = integral(F,0,1);
            HDMR(1) = k;
            return k;
        end
    else
        k = length(M); a = 0;
        for i = 1:k
            a = a + M(i)*(n + 1)^(k - i);
        end
        if HDMR(a) ~= []
            return HDMR(a);
        else
            N = M(1:length(M) - 1)
            k = integral(Kernel(N,F,n),0,1);
            HDMR(a) = k;
            return k;
        end
    end
end

end