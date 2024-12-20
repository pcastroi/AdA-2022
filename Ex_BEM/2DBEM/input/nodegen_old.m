function [xyb,topology,admitt]=nodegen(segments,see);

% [xyb,topology, admitt]=nodegen(segments,see);
%
% Builds geometry arrays from a description of the geometry as
% segments of constant curvature and mesh density. The mesh generated
% consists of quadratic three-node elements. If there are different
% bodies, they are detected and numbered.
%
% Input variables:
%     -segments: array with a segment description per row. The segments
%                must be consecutive in each body and in a clockwise
%                direction. The columns are: x and y of the starting
%                point, x and y of the end point, number of elements in
%                the segment, curvature radius and segment admittance.
%                The sign of the radius indicates whether it bends 
%                outwards (+) or inwards (-). If the radius is 0, 
%                a straight segment is used. 
%  
%     -see:      if is set to 'y' or 'Y', a plot is made with the geometry.
%
% Output variables:
%     -xyb:      node positions, first column is the x-coordinate, second
%                column is y-coordinate, and third column is the body
%                number to which the node belongs to.
%     -topology: each row contains the node numbers (row number in xyb) of the
%                nodes in one element. The last column is the body number the
%                element belongs to.
%     -admitt:   vector admittance for each node

% Susana Quiros y Alpera, Vicente Cutanda 2001

number_seg=size(segments,1);	% number of segments

x=[];
y=[];
bn=[];
xyb=[];
admitt=[];

for m=1:number_seg;	% For every segment
   if m==1;
      body=1;
      first_seg=1;
   else
      if segments(m,1)==segments(m-1,3)&segments(m,2)==segments(m-1,4);
         first_seg=0;
      else
         first_seg=1;
         body=body+1;
         x=x(1:end-1);y=y(1:end-1);bn=bn(1:end-1);admitt=admitt(1:end-1);
      end
   end

   x_start=segments(m,1);
   y_start=segments(m,2);
   x_end=segments(m,3);
   y_end=segments(m,4);
   number_elem=segments(m,5);	% number of elements for the actual segment
   radius=segments(m,6);		% radius=0 for straight segments
   beta=segments(m,7);          % beta== admittance of the segment 
   % if first segment in body
   % jj=1;
   
   number_nod=number_elem*2+first_seg;
   
   % calculate the b's
   bn=[bn body*ones(1,number_nod)];

   % calculate the admittances
   admitt=[admitt ; beta*ones(number_nod,1)];
 
   if radius==0;
      % calculate the x's
      if x_start==x_end;
         x=[x x_start*ones(1,number_nod)];
      else
         step=(x_end-x_start)/(number_elem*2);
         x=[x (x_start+(1-first_seg)*step):step:x_end];
      end
      % calculate the y's
      if y_start==y_end;
         y=[y y_start*ones(1,number_nod)];
   	  else
         step=(y_end-y_start)/(number_elem*2);
         y=[y (y_start+(1-first_seg)*step):step:y_end];
      end
   else % Radius <> 0
      distance=sqrt((y_end-y_start)^2+(x_end-x_start)^2);
      diameter=abs(2*radius);
      if distance>abs(2*radius);	% circle not big enough !!!
         disp(['Wrong input data: Circle not big enough, distance = ' ...
               num2str(distance) ' diameter = ' num2str(diameter)]);
      else
         
         
         phi=angle((x_end-x_start)+i*(y_end-y_start));
         
         x_end2=(x_end-x_start)*cos(phi)+(y_end-y_start)*sin(phi);
         
         y_end2=-(x_end-x_start)*sin(phi)+(y_end-y_start)*cos(phi); % must give zero
         
         radius2=abs(radius)/x_end2;
         
         phi_max=2*asin(0.5/radius2);
         
         step=phi_max/(number_elem*2);
         
         if first_seg
            xi=0;yi=0;
         else
            xi=[];yi=[];
         end
         phi_i=step*[1:(number_elem*2)];
         xi2=[xi 2*radius2*sin(phi_i/2).*cos((phi_max-phi_i)/2)];
         yi2=sign(radius)*[yi 2*radius2*sin(phi_i/2).*sin((phi_max-phi_i)/2)];
         
         xi2=xi2*x_end2;yi2=yi2*x_end2;

         x=[x xi2*cos(phi)-yi2*sin(phi)+x_start];
         y=[y xi2*sin(phi)+yi2*cos(phi)+y_start];
         
         
      end
   end    
   if m==number_seg
      x=x(1:end-1);y=y(1:end-1);bn=bn(1:end-1);admitt=admitt(1:end-1);
   end
end % m loop
% calculate number of elements in each body
xyb=[x' y' bn'];

topology=[];
num_bodies=max(xyb(:,3));
for ii=1:num_bodies
   ind=find(xyb(:,3)==ii);
   topology=[topology ; ind(1:2:end-1) ind(2:2:end) [ind(3:2:end-1);ind(1)] ii*ones(length(ind)/2,1)];
end


if see(1)=='y' | see(1)=='Y'
   figure;
   for bb=1:num_bodies
      eleb=find(topology(:,end)==bb);
      nodb=find(xyb(:,end)==bb);
      plot(xyb(nodb,1),xyb(nodb,2),'r:');
      hold on;
      for nn=1:length(nodb)
         if any(topology(eleb,1)==nodb(nn))|any(topology(eleb,3)==nodb(nn))
            plot(xyb(nodb(nn),1),xyb(nodb(nn),2),'ko');
         else
            plot(xyb(nodb(nn),1),xyb(nodb(nn),2),'k+');
         end
      end
   end
   hold off;
   grid;
   axis equal;
   title(['Nodes = ' num2str(size(xyb,1)) '  Elements = ' num2str(size(topology,1)) ...
         '  Bodies = ' num2str(num_bodies)]);
   xlabel('x');
   ylabel('y');
end


function ang=angle2(a,b);

