function s=count_events_fix_1(T,M)
s=0;
for t=1:T
    if M>t
    s=s+nchoosek(T,t)*nchoosek(M,t)*factorial(t);
    else
        s=s+0;
    end
end

end