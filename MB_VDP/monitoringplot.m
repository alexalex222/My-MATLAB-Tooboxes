clear;
clc;
close all;

testcase1=zeros(600,1);
testcase2=zeros(600,1);


load bcmode1idv1;
testcase1(1:300,1)=classlabel([1:6:900,1001:6:1900],1);
load bcmode3idv11;
testcase1(301:600,1)=classlabel([1:6:900,1001:6:1900],1)+2;

load bcmode1idv14;
testcase2(1:300,1)=classlabel([1:6:900,1001:6:1900],1);
load bcmode3idv13;
testcase2(301:600,1)=classlabel([101:6:1000,1001:6:1900],1)+2;



str1={'\leftarrow -- Base Mode Normal - - \rightarrow'};
str2={'\leftarrow -- Base Mode Faulty - - \rightarrow'};
str3={'\leftarrow -- New Mode Normal - - \rightarrow'};
str4={'\leftarrow -- New Mode Faulty - -\rightarrow'};





figure;
subplot(2,1,1);
plot(testcase1,'o');
hold on;
line([150,150],[0,5],'color','g','linestyle','--');
line([300,300],[0,5],'color','g','linestyle','--');
line([450,450],[0,5],'color','g','linestyle','--');
axis([0 600 0 5]);
set(gca,'YTick',[0:5]);
set(gca,'YTickLabel',{'';'Base Mode Normal';'Base Mode Faulty';'New Mode Normal';'New Mode Fualty'});
title('(a) Test Caes I');
text(75,3.5,str1,'HorizontalAlignment','center','VerticalAlignment','top');
text(225,3.5,str2,'HorizontalAlignment','center','VerticalAlignment','top');
text(375,1.5,str3,'HorizontalAlignment','center','VerticalAlignment','top');
text(525,1.5,str4,'HorizontalAlignment','center','VerticalAlignment','top');
hold off;

subplot(2,1,2);
plot(testcase2,'^');
hold on;
line([150,150],[0,5],'color','g','linestyle','--');
line([300,300],[0,5],'color','g','linestyle','--');
line([450,450],[0,5],'color','g','linestyle','--');
axis([0 600 0 5]);
set(gca,'YTick',[0:5]);
set(gca,'YTickLabel',{'';'Base Mode Normal';'Base Mode Faulty';'New Mode Normal';'New Mode Fualty'});
title('(b) Test Case II');
text(75,3.5,str1,'HorizontalAlignment','center','VerticalAlignment','top');
text(225,3.5,str2,'HorizontalAlignment','center','VerticalAlignment','top');
text(375,1.5,str3,'HorizontalAlignment','center','VerticalAlignment','top');
text(525,1.5,str4,'HorizontalAlignment','center','VerticalAlignment','top');
hold off;




% false alarm rate
falrate1=0;
falrate2=0;
for i=1:150
    if testcase1(i)~=1
        falrate1=falrate1+1;
    end
    
    if testcase2(i)~=1
        falrate2=falrate2+1;
    end
    
end

for i=301:450
    if testcase1(i)~=3
        falrate1=falrate1+1;
    end
    
    if testcase2(i)~=3
        falrate2=falrate2+1;
    end
    

end
falrate1=falrate1/300;
falrate2=falrate2/300;





% fault detection rate
detrate1=300;
detrate2=300;
detrateknn=300;
detratefcm=300;

for i=151:300
    if testcase1(i)==1
    detrate1=detrate1-1;
    end
    
    if testcase2(i)==1
        detrate2=detrate2-1;
    end
    
end

for i=451:600
    if testcase1(i)==3
    detrate1=detrate1-1;
    end
    
    if testcase2(i)==3
        detrate2=detrate2-1;
    end
    
    

end
detrate1=detrate1/300;
detrate2=detrate2/300;
