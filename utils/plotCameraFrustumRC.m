%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2022, Xiao Hu, Francois Lauze, Kim Steenstrup Pedersen
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License. 
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h1 = plotCameraFrustumRC(R,C,K,sz,c,linewidth)
%
    f = K(1,1);
    cu = K(1,3);
    cv = K(2,3);
    if nargin < 4
        sz = 1;
    end
    if nargin < 6
        linewidth = 1;
    end
    
    % compute four rectangle points in camera frame
    z = [sz,sz,sz,sz];
    xx = sz*cu ./ f;
    yy = sz*cv ./ f;
    x = [-xx xx xx -xx];
    y = [-yy -yy yy yy];
    pt_in_cam = [x;y;z];
    camp = C;
    pt_in_world = R'*(pt_in_cam) + repmat(C,1,length(xx));
    % camera center
    plot3(camp(1),camp(2),camp(3),'.','MarkerSize',10,'Color',c);
    % frustum
    for i = 1:4
        if i == 1
            h1 = plot3([camp(1);pt_in_world(1,i)],...
              [camp(2);pt_in_world(2,i)],...
              [camp(3);pt_in_world(3,i)],'-','LineWidth',linewidth,'Color',c);
        else
            plot3([camp(1);pt_in_world(1,i)],...
              [camp(2);pt_in_world(2,i)],...
              [camp(3);pt_in_world(3,i)],'-','LineWidth',linewidth,'Color',c);
        end
    end
    ij = [[1,2];[2,3];[3,4];[4,1]];
    for k = 1:size(ij,1)
        i = ij(k,1);j = ij(k,2);
        plot3([pt_in_world(1,i);pt_in_world(1,j)],...
          [pt_in_world(2,i);pt_in_world(2,j)],...
          [pt_in_world(3,i);pt_in_world(3,j)],'-','LineWidth',linewidth,'Color',c);
    end
end