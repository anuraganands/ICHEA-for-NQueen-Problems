hold on;
X1 = [-100:1:100];
YMatrix1 = -9*X1+6;
YMatrix2 = 9*X1 - 1;
XMatrix3 = 1;
XMatrix4 = 0;

fig1 = gcf; 
axes1 = axes('Parent',fig1); 

ylim(axes1,[-100 100]);
box(axes1,'on');
hold(axes1,'all');
plot(X1,YMatrix1);
plot(X1,YMatrix2,'Parent',axes1);
plot(XMatrix3,X1,'Parent',axes1);
plot(XMatrix4,X1,'Parent',axes1);