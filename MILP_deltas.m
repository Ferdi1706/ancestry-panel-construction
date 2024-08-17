% Implementation of the optimization problem using intlinprog()
function x = MILP_deltas(delta,q,n,m,zweiergruppen,marker_namen,kopplungsgruppen) 
c = delta(:,1);
for i = 2:length(delta(1,:))
    c = c + delta(:,i);
end
c = -(1-q)*c;
c(n+1) = -q;
lb = zeros(n+1,1); ub = ones(n+1,1); ub(n+1) = Inf;
Aeq = ones(1,n+1); Aeq(n+1) = 0; beq = m;
A = ones(length(delta(1,:)),n+1);
for i = 1:length(delta(1,:))
    A(i,1:n) = -delta(:,i);
end
b = zeros(length(delta(1,:)),1);
intcon = 1:n;


% From here on, only construction of constraints for the consideration of coupling
n_rows = length(delta(1,:))+1;

for k = 1:length(kopplungsgruppen)
    gruppe = kopplungsgruppen{k};
    for i = 1:(length(gruppe)-1)
        for j = (i+1):length(gruppe)
            x1 = marker_namen == gruppe(i);
            x2 = marker_namen == gruppe(j);
            A(n_rows, 1:n+1) = zeros(1, n+1);
            A(n_rows, x1) = 1; 
            A(n_rows, x2) = 1; 
            b(n_rows) = 1;
            n_rows = n_rows + 1;
        end
    end
end

k = 1;
while(k<=length(zweiergruppen))
    x1 = marker_namen == zweiergruppen(k);
    x2 = marker_namen == zweiergruppen(k+1);
    A(n_rows,1:n+1) = zeros(1,n+1);
    A(n_rows,x1) = 1; A(n_rows,x2) = 1; b(n_rows) = 1;
    k = k+2;
    n_rows = n_rows + 1;
end

% now find optimal marker combination using intlinprog()
x = intlinprog(c, intcon, A, b, Aeq, beq, lb, ub);
end
