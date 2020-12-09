function v = velocity(xe, k)

   x=xe(:,1);
   gradpx = (1-2*x);
   
   vx = -k*gradpx;
   v = vx;

end