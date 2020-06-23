function p=limitdist(P)
%Obtain the stationary probability distribution
%vector p of an irreducible, recurrent Markov
%chain by state reduction. P is the transition
%probabilities matrix of a discrete-time Markov
%chain or the generator matrix Q.
[ns ms]=size(P);

if P((ns/2+1):ns,1:(ms/2)) == zeros(ns/2, ms/2)
    P = flipud(fliplr(P));
    n=ns;
    while n>1
        n1=n-1;
        s=sum(P(n,1:n1));
        P(1:n1,n)=P(1:n1,n)/s;
        n2=n1;
            while n2>0
                P(1:n1,n2)=P(1:n1,n2)+P(1:n1,n)*P(n,n2);
                n2=n2-1;
            end
        n=n-1;
    end
    %backtracking
    p(1)=1;
    j=2;
    while j<=ns
        j1=j-1;
        p(j)=sum(p(1:j1).*(P(1:j1,j)'));
        j=j+1;
    end
    p=(p/(sum(p)))';

    p = flipud(fliplr(p));
else
    n=ns;
    while n>1
        n1=n-1;
        s=sum(P(n,1:n1));
        P(1:n1,n)=P(1:n1,n)/s;
        n2=n1;
            while n2>0
                P(1:n1,n2)=P(1:n1,n2)+P(1:n1,n)*P(n,n2);
                n2=n2-1;
            end
        n=n-1;
    end
    %backtracking
    p(1)=1;
    j=2;
    while j<=ns
        j1=j-1;
        p(j)=sum(p(1:j1).*(P(1:j1,j)'));
        j=j+1;
    end
    p=(p/(sum(p)))';
end