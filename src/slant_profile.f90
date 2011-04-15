      subroutine slant_profile(Angle_view, Ds_cloud_slant,             &
      Ds_obs_point, DH_bot_intersect, DH_top_intersect)                 
!      given an observation angle, a thickness of the cloud in the      
!      viewing direction (but in the horizontal) and                    
!     a vector of altitudes compute the new vector of altitudes         
!      INPUT=                                                           
!       Ds_cloud_slant (in km)                                          
!       Angle_view in radiants                                          
!      Ds_obs_point in km (positive=out of the clous,                   
!        negative=within the clous)                                     
!       hght_lev=old height vectors                                     
!       hght_lev_new=new height vectors                                 
                                                                        
  use kinds
      implicit none 
      real(kind=dbl) Ds_cloud_slant, Ds_obs_point, angle_view,                 &
      DH_bot_intersect, DH_top_intersect                                
                                                                        
      if (Ds_obs_point.gt.0d0) then 
         DH_bot_intersect = Ds_obs_point * dtan (Angle_view) 
      else 
         DH_bot_intersect = 0.0d0 
      endif 
      DH_top_intersect = (Ds_obs_point + Ds_cloud_slant) * dtan (       &
      Angle_view)                                                       
                                                                        
      return

      end subroutine slant_profile                  
