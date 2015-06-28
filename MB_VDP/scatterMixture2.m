function h = scatterMixture2(data, z)
%function h = scatterMixture(data, z)
    colors = ['g','r','r','r','r','r','r','r','r','r','r'];
    symbols = ['o','x','x','x','x','o','d','v','^','<','>','p','h'];
    for i = min(z):max(z)
        hold on;
        % format = strcat(colors(1+rem(i,numel(colors))),symbols(1+rem(3*i,numel(symbols))));
        format = strcat(colors(i),symbols(i));
        scatter(data(z==i,1),data(z==i,2),format);
    end
    % legend('Base Mode','New Mode','Outlier');
end