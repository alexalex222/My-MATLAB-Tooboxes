
data=T1;

colors = ['r','g','b','c','m','k','y'];
    symbols = ['.','o','x','+','*','s','d','v','^','<','>','p','h'];
%     for i = min(classlabel):max(classlabel)
     for i = min(classlabel):4   
        format = strcat(colors(1+rem(i,numel(colors))),symbols(1+rem(3*i,numel(symbols))));
       
        scatter3(data(classlabel==i,1),data(classlabel==i,2),data(classlabel==i,3),format);
        hold on;
    end

