function [B, div, J] = get_shape_dirv(xi, neldof_u, vtx_coords)
  
      % Calculate the Grad(N) matrix
      if neldof_u == 2        
        GN = [-1 1]/2;
        J = GN*vtx_coords;        % Compute Jacobian
        B = GN/J;                 % compute the derivative of the shape functions
        B1x = B(1,1);
        B2x = B(1,2);
        div = [B1x B2x];
      elseif neldof_u == 3
          
        GN = [xi-1/2 -2*xi xi+1/2];
        J = GN*vtx_coords;        % Compute Jacobian
        B = GN/J;                 % compute the derivative of the shape functions
        B1x = B(1,1);
        B2x = B(1,2);
        B3x = B(1,3);
        div = [B1x B2x B3x];
      else
          error('does not support p>2')
      end


      
end