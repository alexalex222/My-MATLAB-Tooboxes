function h = scatter3dMixture(data, z)

    colors = ['g','b','r','m','c','k','y'];
    symbols = ['.','o','x','+','*','s','d','v','^','<','>','p','h'];
    for i = min(z):max(z)
        
        format = strcat(colors(1+rem(i,numel(colors))),symbols(1+rem(3*i,numel(symbols))));
%         format = strcat(colors(i),symbols(i));
        scatter3(data(z==i,1),data(z==i,2),data(z==i,3),format);
        hold on;
    end
    % legend('Base Mode','New Mode','Outlier');
end