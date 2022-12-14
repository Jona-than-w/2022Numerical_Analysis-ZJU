x=[0.0 0.5 1.0 1.5 2.0 2.5 3.0...
    3.5 4.0 4.5 5.0 5.5 6.0 6.5...
    7.0 7.5 8.0 8.5 9.0 9.5 10.0];
y=[2.9 2.7 4.8 5.3 7.1 7.6 7.7...
 7.6 9.4 9.0 9.6 10.0 10.2 9.7...
 8.3 8.4 9.0 8.3 6.6 6.7 4.1];

A=ones(length(x),3);
for i=1:2
    A(:,i+1)=(x.^i)';
end

[Q,R]=qr(A);
R1=R(1:3,1:3);
c=Q'*y';
c=c(1:3);
a=R1\c;
rn2=norm(R1,2);

f=@(x) a(1)+a(2).*x+a(3).*x.^2;
z=f(x);
plot(x,z);
hold on;
scatter(x,y,'*');