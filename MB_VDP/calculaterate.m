% false alarm rate
falrate=0;
for i=1:1000
    if classlabel(i)~=1
        falrate=falrate+1;
    end

end
falrate=falrate/1000;

% fault detection rate
detrate=1000;
for i=1001:2000
    if classlabel(i)==1
    detrate=detrate-1;
    end
end
detrate=detrate/1000;