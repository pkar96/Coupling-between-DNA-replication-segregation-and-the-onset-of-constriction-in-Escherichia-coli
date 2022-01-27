classdef Cell < matlab.mixin.Copyable
   properties
      t
      v
      rate
      oris
      vNextInit
      vNextDivs
      tNextDivs
      vOfInits
      tOfInits
      oOfInits
      vb
      vi
      vd
      tLastDiv
      tLastInit
      cd_cell
      td_cell
      lb_cell
      ld_cell
      li_cell
      li_t_cell
      rate_cell
      rec_data
   end
   methods
      function obj = Cell(v)
         if nargin == 0
            obj.t = 0;
            obj.v = 1;
            obj.vb = 1;
            obj.vi = 1;
            obj.vd = 1;
            obj.rate = 0;
            obj.vNextInit = 0;
            obj.vNextDivs = [0];
            obj.tNextDivs = [0];
            obj.oris = 2;
            obj.tLastDiv = 0;
            obj.tLastInit = 0;
            obj.vOfInits = [1];
            obj.tOfInits = [0];
            obj.oOfInits = [2];
            obj.td_cell=NaN;
            obj.cd_cell=NaN;
            obj.lb_cell= NaN;
            obj.ld_cell= NaN;
            obj.li_cell= NaN;
            obj.li_t_cell= NaN;
            obj.rate_cell= NaN;
            obj.rec_data=0;
         else
            obj.t = 0;
            obj.v = v;
            obj.vb = v;
            obj.vi = v;
            obj.vd = v;
            obj.rate = 0;
            obj.vNextInit = 0;
            obj.vNextDivs = [0];
            obj.tNextDivs = [0];
            obj.oris = 2;
            obj.tLastDiv = 0;
            obj.tLastInit = 0;
            obj.vOfInits = [v];
            obj.tOfInits = [0];
            obj.oOfInits = [2];
            obj.td_cell= NaN;
            obj.cd_cell= NaN;
            obj.lb_cell= NaN;
            obj.ld_cell= NaN;
            obj.li_cell= NaN;
            obj.li_t_cell= NaN;
            obj.rate_cell= NaN;
            obj.rec_data=0;
         end
      end
   end
end