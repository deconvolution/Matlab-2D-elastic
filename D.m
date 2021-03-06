function u=D(u,direction)
%% first order

switch direction
    case 1
        % backward
        u(2:end,:)=diff(u,1,1);
    case -1
        % forward
        u(1:end-1,:)=diff(u,1,1);
    case 3
        u(:,2:end)=diff(u,1,2);
    case -3
        u(:,1:end-1)=diff(u,1,2);
end

%% third order
%{
switch direction
    case 1
        % backward
        u2=...
            -1/3*ad(u(1:end-3,:),3,-1) ...
            +3/2*ad(u(1:end-2,:),2,-1) ...
            -3*ad(u(1:end-1,:),1,-1) ...
            +11/6*u;
    case -1
        % forward
        u2=-11/6*u ...
            +3*ad(u(2:end,:),1,1) ...
            -3/2*ad(u(3:end,:),2,1) ...
            +1/3*ad(u(4:end,:),3,1);
    case 3
        u2=-1/3*ad(u(:,1:end-3),3,-3) ...
            +3/2*ad(u(:,1:end-2),2,-3) ...
            -3*ad(u(:,1:end-1),1,-3) ...
            +11/6*u;
    case -3
        u2=-11/6*u ...
            +3*ad(u(:,2:end),1,3) ...
            -3/2*ad(u(:,3:end),2,3) ...
            +1/3*ad(u(:,4:end),3,3);
end
%}
%% second order
%{
switch direction
    case 1
        % backward
        u2=1/2*ad(u(1:end-2,:),2,-1) ...
            -2*ad(u(1:end-1,:),1,-1) ...
            +3/2*u;
    case -1
        % forward
        u2=-3/2*u ...
            +2*ad(u(2:end,:),1,1) ...
            -1/2*ad(u(3:end,:),2,1);
    case 3
        u2=1/2*ad(u(:,1:end-2),2,-3) ...
            -2*ad(u(:,1:end-1),1,-3) ...
            +3/2*u;
    case -3
        u2=-3/2*u ...
            +2*ad(u(:,2:end),1,3) ...
            -1/2*ad(u(:,3:end),2,3);
end
%}
%% sixth order
%{
switch direction
    case 1
        % backward
        u2=...
            1/6*ad(u(1:end-6,:),6,-1) ...
            -6/5*ad(u(1:end-5,:),5,-1) ...
            +15/4*ad(u(1:end-4,:),4,-1) ...
            -20/3*ad(u(1:end-3,:),3,-1) ...
            +15/2*ad(u(1:end-2,:),2,-1) ...
            -6*ad(u(1:end-1,:),1,-1) ...
            +49/20*u;
    case -1
        % forward
        u2=-49/20*u ...
            +6*ad(u(2:end,:),1,1) ...
            -15/2*ad(u(3:end,:),2,1) ...
            +20/3*ad(u(4:end,:),3,1) ...
            -15/4*ad(u(5:end,:),4,1) ...
            +6/5*ad(u(6:end,:),5,1) ...
            -1/6*ad(u(7:end,:),6,1);
    case 3
        u2=...
            1/6*ad(u(:,1:end-6),6,-3) ...
            -6/5*ad(u(:,1:end-5),5,-3) ...
            +15/4*ad(u(:,1:end-4),4,-3) ...
            -20/3*ad(u(:,1:end-3),3,-3) ...
            +15/2*ad(u(:,1:end-2),2,-3) ...
            -6*ad(u(:,1:end-1),1,-3) ...
            +49/20*u;
    case -3
        u2=-49/20*u ...
            +6*ad(u(:,2:end),1,3) ...
            -15/2*ad(u(:,3:end),2,3) ...
            +20/3*ad(u(:,4:end),3,3) ...
            -15/4*ad(u(:,5:end),4,3) ...
            +6/5*ad(u(:,6:end),5,3) ...
            -1/6*ad(u(:,7:end),6,3);
end
%}
end

