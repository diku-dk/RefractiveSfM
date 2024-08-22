function Xv = to_virtual_camera_coordinate(X,R,C,Rv,Cv)
    if ~isempty(R)
    	Xc = to_camera_coordinate(X,R,C);
    else
        Xc = X;
    end
    Xv = X;
    for i = 1:length(Rv)
        Xv(:,i) = Rv{i}'*Xc(:,i) - Rv{i}'*Cv{i};
    end
end