module double_moments_module

  use kinds

  implicit none

!Parameter of the drop size distribution when 2moments scheme is used
real (kind=dbl), dimension(4) :: 	gamma_cloud,&
  									gamma_rain, &
  									gamma_ice, &
  									gamma_snow, &
  									gamma_graupel, &
  									gamma_hail

     contains

      subroutine read_double_moments(moments_file)
         character(20) :: dummy
         character(20) :: moments_file

         open(118,file=moments_file)
         !read NU & Mu parameter of the drop size distr. as a function of MASS
         !and Alpha & Beta parameter of Diameter-Mass function
         !gamma_xxx(1)=nu    gamma_xxx(2)=mu     gamma_xxx(3)=alpha      gamma_xxx(4)=beta
         read(118,'(a20,4(x,d13.6))') dummy,gamma_cloud
         read(118,'(a20,4(x,d13.6))') dummy,gamma_rain
         read(118,'(a20,4(x,d13.6))') dummy,gamma_ice
         read(118,'(a20,4(x,d13.6))') dummy,gamma_snow
         read(118,'(a20,4(x,d13.6))') dummy,gamma_graupel
         read(118,'(a20,4(x,d13.6))') dummy,gamma_hail
         close(118)

         return
      end subroutine read_double_moments

end module double_moments_module
