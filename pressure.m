function p = pressure(xe)

   x=xe(:,1);
   
   p = x.*(1-x);

end