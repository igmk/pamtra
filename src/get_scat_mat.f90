subroutine get_scat_mat(layer,nstokes, nummu,&
scatter_matrix,&
extinct_matrix, emis_vector)
    use kinds
    use vars_atmosphere, only: rt4scatter_matrix, rt4ext_matrix, rt4emis_vec

    implicit none

    integer :: layer, nstokes, nummu, j, l, i, i1, i2, j1, j2, l1, l2
    real*8  scatter_matrix(nstokes,nummu,nstokes,nummu,4)
    real*8  extinct_matrix(nstokes,nstokes,nummu,2)
    real*8  emis_vector(nstokes,nummu,2)

    do l1 = 1, 2
        do j1 = 1, nummu
            do l2 = 1, 2
                l = 2*(l2-1)+l1
                do j2 = 1, nummu
                    do i2 = 1, nstokes
                        do i1=1,nstokes
                            scatter_matrix(i2,j2,i1,j1,l) = rt4scatter_matrix(layer,i2,j2,i1,j1,l)
                        end do
                    enddo

                enddo
            enddo
        enddo
    enddo

    do l = 1, 2
        do j = 1, nummu
            do i2 = 1, nstokes
                do i1 = 1, nstokes
                    extinct_matrix(i2,i1,j,l) = rt4ext_matrix(layer,i2,i1,j,l)
                end do
            enddo
        enddo
    enddo

    do l = 1, 2
        do j = 1, nummu
            do i = 1,nstokes
                emis_vector(i,j,l) = rt4emis_vec(layer,i,j,l)
            end do
        enddo
    enddo


    return
end subroutine get_scat_mat
