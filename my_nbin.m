function output = my_nbin(min,max,M,k,cdf)
if cdf==true
    %For cdf = true, gives cdf value for x between min and max, for each
    %element of a vector M
    output = zeros(1,length(M)); %Row vector
    for i = 1:length(M)
        X = [min:max];
        output(i) = sum(gamma(k+X)./(gamma(X+1).*gamma(k)).*(M(i)./(M(i)+k)).^(X).*(1+(M(i)./k)).^(-k));
    end
end
if cdf==false
    %If cdf = false, outputs pdf value for each x in range min to max, for
    %1 vlaue of M.
    X = [min:max];
    output = gamma(k+X)./(gamma(X+1).*gamma(k)).*(M./(M+k)).^(X).*(1+(M./k)).^(-k);
end