function d2q=momenluc(q,dq,to,fd)
d1=0.1397;
a2=0.4318;
d3=0.4337;
d5=0.05588;
q1=q(1);
q2=q(2);
q3=q(3);
q4=q(4);
q5=q(5);
q6=q(6);
A1=[cos(q1) 0 sin(q1) d1*-sin(q1);
    sin(q1) 0 -cos(q1) d1*cos(q1);
    0 1 0 0;
    0 0 0 1];
A2=[cos(q2) -sin(q2) 0 a2*cos(q2);
    sin(q2) cos(q2) 0 a2*sin(q2);
    0 0 1 0;
    0 0 0 1];
A3=[cos(q3) 0 -sin(q3) d3*-sin(q3);
    sin(q3) 0 cos(q3) d3*cos(q3);
    0 -1 0 0;
    0 0 0 1];
A4=[cos(q4) 0 sin(q4) 0;
    sin(q4) 0 -cos(q4) 0;
    0 1 0 0;
    0 0 0 1];
A5=[cos(q5) 0 -sin(q5) -d5*sin(q5);
    sin(q5) 0 cos(q5) d5*cos(q5);
    0 -1 0 0;
    0 0 0 1];
A6=[cos(q6) -sin(q6) 0 0;
    sin(q6) cos(q6) 0 0;
    0 0 1 0;
    0 0 0 1];
A=zeros(4,4,6);
A(:,:,1)=A1;
A(:,:,2)=A2;
A(:,:,3)=A3;
A(:,:,4)=A4;
A(:,:,5)=A5;
A(:,:,6)=A6;
m1=17.4;
m2=4.8;
m3=0.82;
m4=0.34;
m5=0.09;
m6=0.05;
J=zeros(4,4,6);
J(:,:,1)=[0 0 0 0;
    0 0 0 0;
    0 0 m1*d1^2/3 m1*d1/2;
    0 0 m1*d1/2 m1];
J(:,:,2)=[m2*a2^2/3 0 0 -m2*a2/2;
    0 0 0 0;
    0 0 0 0;
    -m2*a2/2 0 0 m2];
J(:,:,3)=[0 0 0 0;
    0 0 0 0;
    0 0 m3*d3^2/3 -m3*d3/2;
    0 0 -m3*d3/2 m3];
J(:,:,4)=[m4*d5^2/4 0 0 0;
    0 m4*d5^2/4 0 0;
    0 0 0 0;
    0 0 0 m4];
J(:,:,5)=[m5*d5^2/4 0 0 0;
    0 0 0 0;
    0 0 5*m5*d5^2/4 -m5*d5;
    0 0 -m5*d5 m5];

J(:,:,6)=[m6*0.05^2 0 0 m6*0.05;
    0 0 0 0;
    0 0 0 0;
    m6*0.05 0 0 m6];
U=zeros(4,4,6,6);
g=[0 0 -9.81 1]';
r=[0 -a2/2 0 0 0 0.05;
    0 0 0 0 0 0;
    d1/2 0 -d3/2 0 -d5 0;
    1 1 1 1 1 1];
for i=1:6
    for j=1:i
        U(:,:,i,j)=eye(4);
        for k=1:j-1
            U(:,:,i,j)=U(:,:,i,j)*A(:,:,k);
        end
        U(:,:,i,j)=U(:,:,i,j)*[0 -1 0 0;1 0 0 0;0 0 0 0;0 0 0 0];
        for k=j:i
            U(:,:,i,j)=U(:,:,i,j)*A(:,:,k);
        end
    end
end
M=zeros(6);
for i=1:6
    for j=i:6
        for k=max(i,j):6
            M(i,j)=M(i,j)+trace(U(:,:,k,j)*J(:,:,k)*U(:,:,k,i)');
        end
        M(j,i)=M(i,j);
    end
end
G=zeros(6,1);
m=[m1;m2;m3;m4;m5;m6];
for i=1:6
    for j=i:6
        G(i)=G(i)+(-m(j)*g'*U(:,:,j,i)*r(:,i));
    end
end
U1=zeros(4,4,6,6,6);
for i=1:6
    for j=1:i
        for k=1:j
            U1(:,:,i,j,k)=eye(4);
            for k1=1:k-1
                U1(:,:,i,j,k)=U1(:,:,i,j,k)*A(:,:,i);
            end
            U1(:,:,i,j,k)=U1(:,:,i,j,k)*[0 -1 0 0;1 0 0 0;0 0 0 0;0 0 0 0];
            for k1=k:j-1
                U1(:,:,i,j,k)=U1(:,:,i,j,k)*A(:,:,i);
            end
            U1(:,:,i,j,k)=U1(:,:,i,j,k)*[0 -1 0 0;1 0 0 0;0 0 0 0;0 0 0 0];
            for k1=j:i
                U1(:,:,i,j,k)=U1(:,:,i,j,k)*A(:,:,i);
            end
            U1(:,:,i,k,j)=U1(:,:,i,j,k);
        end
    end
end
h=zeros(6,6,6);
for i=1:6
    for j=1:6
        for k=1:6
            for k1=max([i,j,k]):6
                h(j,k,i)=h(j,k,i)+trace(U1(:,:,k1,j,k)*J(k1)*U(:,:,k1,i));
            end
        end
    end
end
V=zeros(6,1);
for i=1:6
    V(i)=dq'*h(:,:,i)*dq;
end
d2q=M\(to-V-G-fd);