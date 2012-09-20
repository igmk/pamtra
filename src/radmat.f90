!      This file contains a number of simple matrix manipulation routine
!      All of the routines operate on REAL*8 matrices.                  
!      The size is passed as the first parameters:  N rows by M columns.
!      If there is no M parameter then a N by N square matrix is assumed
!                                                                       
!    CALL MCOPY (N,M, X, Y)        Copies matrix: Y = X                 
!                                                                       
!    CALL MADD (N,M, X, Y, Z)      Adds matrices: Z = X + Y             
!                                    Z can be X or Y                    
!    CALL MSUB (N,M, X, Y, Z)      Subtracts matrices: Z = X - Y        
!                                    Z can be X or Y                    
!    CALL MSCALARMULT (N,M, C, X, Y) Scalar multiply: Y = C*X           
!                                    C is real scalar; Y can be X       
!    CALL MZERO (N,M, X)           Zeros all elements in X              
!                                                                       
!    CALL MDIAG (N, V, X)          Formats a vector (V) into a diagonal 
!                                    matrix (X)                         
!    CALL MIDENTITY (N, X)         Make identity matrix in X            
!                                                   t                   
!    CALL MTRANSPOSE (N,M, X, Y)   Transposes: Y = X ;   Y cannot be X  
!                                                                       
!    CALL MMULT (N,M,L, X, Y, Z)   Matrix multiply: Z = X*Y             
!                                    X is N by M and Y is M by L.       
!                                    Z cannot be X or Y                 
!                                                          -1           
!    CALL MINVERT (N, X, Y)        Matrix inversion:  Y = X             
!				     X gets LU decomposition.                                      
      SUBROUTINE MCOPY (N, M, MATRIX1, MATRIX2) 
      use kinds
      INTEGER N, M, I 
      REAL(kind=dbl) MATRIX1 (1), MATRIX2 (1) 
                                                                        
      DO 100 I = 1, N * M 
         MATRIX2 (I) = MATRIX1 (I) 
  100 END DO 
                                                                        
      RETURN 
      END SUBROUTINE MCOPY                          
                                                                        
                                                                        
      SUBROUTINE MADD (N, M, MATRIX1, MATRIX2, MATRIX3) 
      use kinds
      INTEGER N, M, I 
      REAL(kind=dbl) MATRIX1 (1), MATRIX2 (1), MATRIX3 (1) 
                                                                        
      DO 100 I = 1, N * M 
         MATRIX3 (I) = MATRIX1 (I) + MATRIX2 (I) 
  100 END DO 
                                                                        
      RETURN 
      END SUBROUTINE MADD                           
                                                                        
                                                                        
      SUBROUTINE MSUB (N, M, MATRIX1, MATRIX2, MATRIX3) 
      use kinds
      INTEGER N, M, I 
      REAL(kind=dbl) MATRIX1 (1), MATRIX2 (1), MATRIX3 (1) 
                                                                        
      DO 100 I = 1, N * M 
         MATRIX3 (I) = MATRIX1 (I) - MATRIX2 (I) 
  100 END DO 
                                                                        
      RETURN 
      END SUBROUTINE MSUB                           
                                                                        
                                                                        
      SUBROUTINE MSCALARMULT (N, M, C, MATRIX1, MATRIX2) 
      use kinds
      INTEGER N, M, I 
      REAL(kind=dbl) C, MATRIX1 (1), MATRIX2 (1) 
                                                                        
      DO 100 I = 1, N * M 
         MATRIX2 (I) = C * MATRIX1 (I) 
  100 END DO 
                                                                        
      RETURN 
      END SUBROUTINE MSCALARMULT                    
                                                                        
                                                                        
      SUBROUTINE MZERO (N, M, MATRIX1) 
      use kinds
      INTEGER N, M, I 
      REAL(kind=dbl) MATRIX1 (1) 
                                                                        
      DO 100 I = 1, N * M 
         MATRIX1 (I) = 0.0D0 
  100 END DO 
                                                                        
      RETURN 
      END SUBROUTINE MZERO                          
                                                                        
                                                                        
      SUBROUTINE MDIAG (N, VECTOR, MATRIX) 
      use kinds
      INTEGER N, I, J 
      REAL(kind=dbl) VECTOR (1), MATRIX (N, N) 
                                                                        
      DO 110 I = 1, N 
         DO 100 J = 1, N 
            MATRIX (I, J) = 0.0 
  100    END DO 
         MATRIX (I, I) = VECTOR (I) 
  110 END DO 
                                                                        
      RETURN 
      END SUBROUTINE MDIAG                          
                                                                        
                                                                        
      SUBROUTINE MIDENTITY (N, MATRIX) 
      use kinds
      INTEGER N, I, J 
      REAL(kind=dbl) MATRIX (N, N) 
                                                                        
      DO 110 I = 1, N 
         DO 100 J = 1, N 
            MATRIX (I, J) = 0.0 
  100    END DO 
         MATRIX (I, I) = 1.0 
  110 END DO 
                                                                        
      RETURN 
      END SUBROUTINE MIDENTITY                      
                                                                        
                                                                        
      SUBROUTINE MTRANSPOSE (N, M, MATRIX1, MATRIX2) 
      use kinds
      INTEGER N, M, I, J 
      REAL(kind=dbl) MATRIX1 (N, M), MATRIX2 (M, N) 
                                                                        
      DO 100 I = 1, N 
         DO 100 J = 1, M 
            MATRIX2 (I, J) = MATRIX1 (J, I) 
  100 CONTINUE 
                                                                        
      RETURN 
      END SUBROUTINE MTRANSPOSE                     
                                                                        
                                                                        
      SUBROUTINE MMULT (N, M, L, MATRIX1, MATRIX2, MATRIX3) 
      use kinds
      INTEGER N, M, L, I, J, K 
      REAL(kind=dbl) MATRIX1 (N, M), MATRIX2 (M, L), MATRIX3 (N, L), SUM 
                                                                        
      DO 200 I = 1, N 
         DO 200 J = 1, L 
            SUM = 0.0 
            DO 100 K = 1, M 
               SUM = SUM + MATRIX1 (I, K) * MATRIX2 (K, J) 
  100       END DO 
            MATRIX3 (I, J) = SUM 
  200 CONTINUE 
      RETURN 
      END SUBROUTINE MMULT                          
                                                                        
                                                                        
!       SUBROUTINE MINVERT_OLD (N, MATRIX1, MATRIX2) 
!       use kinds
!       INTEGER N 
!       REAL(kind=dbl) MATRIX1 (N, N), MATRIX2 (N, N) 
!       INTEGER NMAX 
!       PARAMETER (NMAX = 256) 
!       INTEGER I, J, INDX (NMAX), IZ 
!       REAL(kind=dbl) DET (2), WORK (NMAX) 
!                                                                         
!       IF (N.GT.NMAX) THEN 
!       WRITE ( * , '(1X,A,I3)') 'Exceeded maximum matrix size for inversi&
!      &on.  Max = ', NMAX                                                
!          STOP 
!       ENDIF 
!       CALL DGEFA (MATRIX1, N, N, INDX, IZ) 
!       IF (IZ.GT.0) THEN 
!       WRITE ( * , '(1X,A,I3)') 'Encountered a zero pivot at element ', I&
!      &Z                                                                 
!          STOP 
!       ENDIF 
!       DO 100 I = 1, N 
!          DO 110 J = 1, N 
!             MATRIX2 (I, J) = MATRIX1 (I, J) 
!   110    END DO 
!   100 END DO 
!       CALL DGEDI (MATRIX2, N, N, INDX, DET, WORK, 1) 
!                                                                         
!       RETURN 
!       END SUBROUTINE MINVERT_OLD
                                                                        

      SUBROUTINE MINVERT (N, MATRIX1, MATRIX2) 
      use kinds
      INTEGER N 
      REAL(kind=dbl) MATRIX1 (N, N), MATRIX2 (N, N) 
      INTEGER IPIV (N), WORK (N*N)
      INTEGER INFO, LDA, M, LWORK
      ! method found at
      ! http://vibrationdata.com/python-wiki/index.php?title=Matrix_Inversion
      external DGETRF
      external DGETRI                                                                      

!     DGETRF computes an LU factorization of a general M-by-N matrix A
!     using partial pivoting with row interchanges.

      LWORK = N*N
      M=N
      LDA=N

!      Store MATRIX1 in MATRIX2 to prevent it from being overwritten by LAPACK
      MATRIX2 = MATRIX1
      CALL DGETRF( M, N, MATRIX2, LDA, IPIV, INFO )
      IF(INFO.LT.0)THEN
	  PRINT '(" LU decomposition:  illegal value ")'
	  STOP
	ENDIF

      CALL DGETRI(N, MATRIX2, N, IPIV, WORK, LWORK, INFO)
      IF (info.NE.0) THEN
         stop 'Matrix inversion failed!'
      ENDIF

      RETURN
      END SUBROUTINE MINVERT                




                                                                        
      SUBROUTINE dgefa (a, lda, n, ipvt, info) 
  use kinds
      INTEGER lda, n, ipvt (1), info 
      real(kind=dbl) a (lda, 1) 
!                                                                       
!     dgefa factors a double precision matrix by gaussian elimination.  
!                                                                       
!     dgefa is usually called by dgeco, but it can be called            
!     directly with a saving in time if  rcond  is not needed.          
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .                   
!                                                                       
!     on entry                                                          
!                                                                       
!        a       double precision(lda, n)                               
!                the matrix to be factored.                             
!                                                                       
!        lda     integer                                                
!                the leading dimension of the array  a .                
!                                                                       
!        n       integer                                                
!                the order of the matrix  a .                           
!                                                                       
!     on return                                                         
!                                                                       
!        a       an upper triangular matrix and the multipliers         
!                which were used to obtain it.                          
!                the factorization can be written  a = l*u  where       
!                l  is a product of permutation and unit lower          
!                triangular matrices and  u  is upper triangular.       
!                                                                       
!        ipvt    integer(n)                                             
!                an integer vector of pivot indices.                    
!                                                                       
!        info    integer                                                
!                = 0  normal value.                                     
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error       
!                     condition for this subroutine, but it does        
!                     indicate that dgesl or dgedi will divide by zero  
!                     if called.  use  rcond  in dgeco for a reliable   
!                     indication of singularity.                        
!                                                                       
!     linpack. this version dated 08/14/78 .                            
!     cleve moler, university of new mexico, argonne national lab.      
!                                                                       
!     subroutines and functions                                         
!                                                                       
!     blas daxpy,dscal,idamax                                           
!                                                                       
!     internal variables                                                
!                   
      real(kind=dbl) t 
      INTEGER idamax, j, k, kp1, l, nm1 
!                                                                       
!                                                                       
!     gaussian elimination with partial pivoting                        
!                                                                       
      info = 0 
      nm1 = n - 1 
      IF (nm1.lt.1) goto 70 
      DO 60 k = 1, nm1 
         kp1 = k + 1 
!                                                                       
!        find l = pivot index                                           
!                                                                       
         l = idamax (n - k + 1, a (k, k), 1) + k - 1 
         ipvt (k) = l 
!                                                                       
!        zero pivot implies this column already triangularized          
!                                                                       
         IF (a (l, k) .eq.0.0d0) goto 40 
!                                                                       
!           interchange if necessary                                    
!                                                                       
         IF (l.eq.k) goto 10 
         t = a (l, k) 
         a (l, k) = a (k, k) 
         a (k, k) = t 
   10    CONTINUE 
!                                                                       
!           compute multipliers                                         
!                                                                       
         t = - 1.0d0 / a (k, k) 
         CALL dscal (n - k, t, a (k + 1, k), 1) 
!                                                                       
!           row elimination with column indexing                        
!                                                                       
         DO 30 j = kp1, n 
            t = a (l, j) 
            IF (l.eq.k) goto 20 
            a (l, j) = a (k, j) 
            a (k, j) = t 
   20       CONTINUE 
            CALL daxpy (n - k, t, a (k + 1, k), 1, a (k + 1, j),        &
            1)                                                          
   30    END DO 
         GOTO 50 
   40    CONTINUE 
         info = k 
   50    CONTINUE 
   60 END DO 
   70 CONTINUE 
      ipvt (n) = n 
      IF (a (n, n) .eq.0.0d0) info = n 
      RETURN 
      END SUBROUTINE dgefa                          
                                                                        
                                                                        
      SUBROUTINE dgedi (a, lda, n, ipvt, det, work, job) 
  use kinds
      INTEGER lda, n, ipvt (1), job 
      real(kind=dbl) a (lda, 1), det (2), work (1) 
!                                                                       
!     dgedi computes the determinant and inverse of a matrix            
!     using the factors computed by dgeco or dgefa.                     
!                                                                       
!     on entry                                                          
!                                                                       
!        a       double precision(lda, n)                               
!                the output from dgeco or dgefa.                        
!                                                                       
!        lda     integer                                                
!                the leading dimension of the array  a .                
!                                                                       
!        n       integer                                                
!                the order of the matrix  a .                           
!                                                                       
!        ipvt    integer(n)                                             
!                the pivot vector from dgeco or dgefa.                  
!                                                                       
!        work    double precision(n)                                    
!                work vector.  contents destroyed.                      
!                                                                       
!        job     integer                                                
!                = 11   both determinant and inverse.                   
!                = 01   inverse only.                                   
!                = 10   determinant only.                               
!                                                                       
!     on return                                                         
!                                                                       
!        a       inverse of original matrix if requested.               
!                otherwise unchanged.                                   
!                                                                       
!        det     double precision(2)                                    
!                determinant of original matrix if requested.           
!                otherwise not referenced.                              
!                determinant = det(1) * 10.0**det(2)                    
!                with  1.0 .le. dabs(det(1)) .lt. 10.0                  
!                or  det(1) .eq. 0.0 .                                  
!                                                                       
!     error condition                                                   
!                                                                       
!        a division by zero will occur if the input factor contains     
!        a zero on the diagonal and the inverse is requested.           
!        it will not occur if the subroutines are called correctly      
!        and if dgeco has set rcond .gt. 0.0 or dgefa has set           
!        info .eq. 0 .                                                  
!                                                                       
!     linpack. this version dated 08/14/78 .                            
!     cleve moler, university of new mexico, argonne national lab.      
!                                                                       
!     subroutines and functions                                         
!                                                                       
!     blas daxpy,dscal,dswap                                            
!     fortran dabs,mod                                                  
!                                                                       
!     internal variables                                                
!                                                                       
      real(kind=dbl) t 
      real(kind=dbl) ten 
      INTEGER i, j, k, kb, kp1, l, nm1 
!                                                                       
!                                                                       
!     compute determinant                                               
!                                                                       
      IF (job / 10.eq.0) goto 70 
      det (1) = 1.0d0 
      det (2) = 0.0d0 
      ten = 10.0d0 
      DO 50 i = 1, n 
         IF (ipvt (i) .ne.i) det (1) = - det (1) 
         det (1) = a (i, i) * det (1) 
!        ...exit                                                        
         IF (det (1) .eq.0.0d0) goto 60 
   10    IF (dabs (det (1) ) .ge.1.0d0) goto 20 
         det (1) = ten * det (1) 
         det (2) = det (2) - 1.0d0 
         GOTO 10 
   20    CONTINUE 
   30    IF (dabs (det (1) ) .lt.ten) goto 40 
         det (1) = det (1) / ten 
         det (2) = det (2) + 1.0d0 
         GOTO 30 
   40    CONTINUE 
   50 END DO 
   60 CONTINUE 
   70 CONTINUE 
!                                                                       
!     compute inverse(u)                                                
!                                                                       
      IF (mod (job, 10) .eq.0) goto 150 
      DO 100 k = 1, n 
         a (k, k) = 1.0d0 / a (k, k) 
         t = - a (k, k) 
         CALL dscal (k - 1, t, a (1, k), 1) 
         kp1 = k + 1 
         IF (n.lt.kp1) goto 90 
         DO 80 j = kp1, n 
            t = a (k, j) 
            a (k, j) = 0.0d0 
            CALL daxpy (k, t, a (1, k), 1, a (1, j), 1) 
   80    END DO 
   90    CONTINUE 
  100 END DO 
!                                                                       
!        form inverse(u)*inverse(l)                                     
!                                                                       
      nm1 = n - 1 
      IF (nm1.lt.1) goto 140 
      DO 130 kb = 1, nm1 
         k = n - kb 
         kp1 = k + 1 
         DO 110 i = kp1, n 
            work (i) = a (i, k) 
            a (i, k) = 0.0d0 
  110    END DO 
         DO 120 j = kp1, n 
            t = work (j) 
            CALL daxpy (n, t, a (1, j), 1, a (1, k), 1) 
  120    END DO 
         l = ipvt (k) 
         IF (l.ne.k) call dswap (n, a (1, k), 1, a (1, l), 1) 
  130 END DO 
  140 CONTINUE 
  150 CONTINUE 
      RETURN 
      END SUBROUTINE dgedi                          
                                                                        
                                                                        
      SUBROUTINE daxpy (n, da, dx, incx, dy, incy) 
!                                                                       
!     constant times a vector plus a vector.                            
!     uses unrolled loops for increments equal to one.                  
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
  use kinds
      real(kind=dbl)  dx (1), dy (1), da 
      INTEGER i, incx, incy, ix, iy, m, mp1, n 
!                                                                       
      IF (n.le.0) return 
      IF (da.eq.0.0d0) return 
      IF (incx.eq.1.and.incy.eq.1) goto 20 
!                                                                       
!        code for unequal increments or equal increments                
!          not equal to 1                                               
!                                                                       
      ix = 1 
      iy = 1 
      IF (incx.lt.0) ix = ( - n + 1) * incx + 1 
      IF (incy.lt.0) iy = ( - n + 1) * incy + 1 
      DO 10 i = 1, n 
         dy (iy) = dy (iy) + da * dx (ix) 
         ix = ix + incx 
         iy = iy + incy 
   10 END DO 
      RETURN 
!                                                                       
!        code for both increments equal to 1                            
!                                                                       
!                                                                       
!        clean-up loop                                                  
!                                                                       
   20 m = mod (n, 4) 
      IF (m.eq.0) goto 40 
      DO 30 i = 1, m 
         dy (i) = dy (i) + da * dx (i) 
   30 END DO 
      IF (n.lt.4) return 
   40 mp1 = m + 1 
      DO 50 i = mp1, n, 4 
         dy (i) = dy (i) + da * dx (i) 
         dy (i + 1) = dy (i + 1) + da * dx (i + 1) 
         dy (i + 2) = dy (i + 2) + da * dx (i + 2) 
         dy (i + 3) = dy (i + 3) + da * dx (i + 3) 
   50 END DO 
      RETURN 
      END SUBROUTINE daxpy                          
                                                                        
                                                                        
      SUBROUTINE dscal (n, da, dx, incx) 
!                                                                       
!     scales a vector by a constant.                                    
!     uses unrolled loops for increment equal to one.                   
!     jack dongarra, linpack, 3/11/78.                                  
!     modified 3/93 to return if incx .le. 0.                           
!                                                                       
  use kinds
      real(kind=dbl) da, dx (1) 
      INTEGER i, incx, m, mp1, n, nincx 
!                                                                       
      IF (n.le.0.or.incx.le.0) return 
      IF (incx.eq.1) goto 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      nincx = n * incx 
      DO 10 i = 1, nincx, incx 
         dx (i) = da * dx (i) 
   10 END DO 
      RETURN 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
!                                                                       
!        clean-up loop                                                  
!                                                                       
   20 m = mod (n, 5) 
      IF (m.eq.0) goto 40 
      DO 30 i = 1, m 
         dx (i) = da * dx (i) 
   30 END DO 
      IF (n.lt.5) return 
   40 mp1 = m + 1 
      DO 50 i = mp1, n, 5 
         dx (i) = da * dx (i) 
         dx (i + 1) = da * dx (i + 1) 
         dx (i + 2) = da * dx (i + 2) 
         dx (i + 3) = da * dx (i + 3) 
         dx (i + 4) = da * dx (i + 4) 
   50 END DO 
      RETURN 
      END SUBROUTINE dscal                          
                                                                        
                                                                        
      INTEGER function idamax (n, dx, incx) 
!                                                                       
!     finds the index of element having max. absolute value.            
!     jack dongarra, linpack, 3/11/78.                                  
!     modified 3/93 to return if incx .le. 0.                           
!                                                                       
  use kinds
      real(kind=dbl) :: dx (1), dmax 
      INTEGER i, incx, ix, n 
!                                                                       
      idamax = 0 
      IF (n.lt.1.or.incx.le.0) return 
      idamax = 1 
      IF (n.eq.1) return 
      IF (incx.eq.1) goto 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      ix = 1 
      dmax = dabs (dx (1) ) 
      ix = ix + incx 
      DO 10 i = 2, n 
         IF (dabs (dx (ix) ) .le.dmax) goto 5 
         idamax = i 
         dmax = dabs (dx (ix) ) 
    5    ix = ix + incx 
   10 END DO 
      RETURN 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
   20 dmax = dabs (dx (1) ) 
      DO i = 2, n 
         IF (dabs (dx (i) ) .le. dmax) exit
         idamax = i 
         dmax = dabs (dx (i) ) 
      END DO
      RETURN 
      END FUNCTION idamax                           
                                                                        
                                                                        
      SUBROUTINE dswap (n, dx, incx, dy, incy) 
!                                                                       
!     interchanges two vectors.                                         
!     uses unrolled loops for increments equal one.                     
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
  use kinds
      real(kind=dbl) :: dx (1), dy (1), dtemp 
      INTEGER i, incx, incy, ix, iy, m, mp1, n 
!                                                                       
      IF (n.le.0) return 
      IF (incx.eq.1.and.incy.eq.1) goto 20 
!                                                                       
!       code for unequal increments or equal increments not equal       
!         to 1                                                          
!                                                                       
      ix = 1 
      iy = 1 
      IF (incx.lt.0) ix = ( - n + 1) * incx + 1 
      IF (incy.lt.0) iy = ( - n + 1) * incy + 1 
      DO 10 i = 1, n 
         dtemp = dx (ix) 
         dx (ix) = dy (iy) 
         dy (iy) = dtemp 
         ix = ix + incx 
         iy = iy + incy 
   10 END DO 
      RETURN 
!                                                                       
!       code for both increments equal to 1                             
!                                                                       
!                                                                       
!       clean-up loop                                                   
!                                                                       
   20 m = mod (n, 3) 
      IF (m.eq.0) goto 40 
      DO 30 i = 1, m 
         dtemp = dx (i) 
         dx (i) = dy (i) 
         dy (i) = dtemp 
   30 END DO 
      IF (n.lt.3) return 
   40 mp1 = m + 1 
      DO 50 i = mp1, n, 3 
         dtemp = dx (i) 
         dx (i) = dy (i) 
         dy (i) = dtemp 
         dtemp = dx (i + 1) 
         dx (i + 1) = dy (i + 1) 
         dy (i + 1) = dtemp 
         dtemp = dx (i + 2) 
         dx (i + 2) = dy (i + 2) 
         dy (i + 2) = dtemp 
   50 END DO 
      RETURN 
      END SUBROUTINE dswap                          
