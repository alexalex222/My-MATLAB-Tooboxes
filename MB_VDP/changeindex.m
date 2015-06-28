for i=1:1000
    if classlabel(i)>2
        classlabel(i)=2;
    end
end
for i=1:2000
    if classlabel(i)==2;
       classlabel(i)=1;
    else classlabel(i)=2;
    end
end