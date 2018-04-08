function id = roullete(S)
n = length(S);
S = S/sum(S);
[s,count] = sort(S,'descend');
for i = 1:n
    if i == 1
        tmp(i) = s(i);
    else
        tmp(i) = tmp(i-1) + s(i);
    end
end
tmp = tmp.^(1/2);
r = rand;
for i = 1:n
    if r < tmp(i)
        id = i;
        break;
    end
end
id = count(id);
% plot([1:n],tmp,'-b');
end
