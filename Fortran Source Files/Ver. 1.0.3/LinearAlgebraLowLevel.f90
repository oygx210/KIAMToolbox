       module LinearAlgebraLowLevel

    contains

       ! About function, open for F2PY wrapping
       function linearalgebralowlevel_about() result(i)
            integer :: i
            i = 0
            write(*,*) 'LinearAlgebraLowLevel: Version 1.0.'
            write(*,*) 'Author: Maksim Shirobokov.'
            write(*,*) 'Date: 19.01.2022.'
       end function linearalgebralowlevel_about

       SUBROUTINE dgetrs( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

       CHARACTER          TRANS
       INTEGER            INFO, LDA, LDB, N, NRHS

       INTEGER            IPIV( * )
       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )

       DOUBLE PRECISION   ONE
       parameter( one = 1.0d+0 )

       LOGICAL            NOTRAN

       !LOGICAL            LSAME
       !EXTERNAL           lsame

       INTRINSIC          max

       info = 0
       notran = lsame( trans, 'N' )
       IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) .AND. .NOT. lsame( trans, 'C' ) ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( nrhs.LT.0 ) THEN
          info = -3
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -5
       ELSE IF( ldb.LT.max( 1, n ) ) THEN
          info = -8
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'DGETRS', -info )
          RETURN
       END IF

       IF( n.EQ.0 .OR. nrhs.EQ.0 ) RETURN

       IF( notran ) THEN

          CALL dlaswp( nrhs, b, ldb, 1, n, ipiv, 1 )

          CALL dtrsm( 'Left', 'Lower', 'No transpose', 'Unit', n, nrhs, one, a, lda, b, ldb )

          CALL dtrsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n, nrhs, one, a, lda, b, ldb )
       ELSE

          CALL dtrsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs, one, a, lda, b, ldb )

          CALL dtrsm( 'Left', 'Lower', 'Transpose', 'Unit', n, nrhs, one, a, lda, b, ldb )

          CALL dlaswp( nrhs, b, ldb, 1, n, ipiv, -1 )
       END IF

       RETURN

       END SUBROUTINE dgetrs

       SUBROUTINE dgetrf ( M, N, A, LDA, IPIV, INFO)

       INTEGER            INFO, LDA, M, N

       INTEGER            IPIV( * )
       DOUBLE PRECISION   A( LDA, * )

       DOUBLE PRECISION   ONE
       parameter( one = 1.0d+0 )

       INTEGER            I, IINFO, J, JB, NB

       !INTEGER            ILAENV

       INTRINSIC          max, min

       info = 0
       IF( m.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = -4
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'DGETRF', -info )
          RETURN
       END IF

       IF( m.EQ.0 .OR. n.EQ.0 ) RETURN

       nb = ilaenv( 1, 'DGETRF', ' ', m, n, -1, -1 )
       IF( nb.LE.1 .OR. nb.GE.min( m, n ) ) THEN

          CALL dgetf2( m, n, a, lda, ipiv, info )
       ELSE

          DO 20 j = 1, min( m, n ), nb
             jb = min( min( m, n )-j+1, nb )

             CALL dgemm( 'No transpose', 'No transpose', m-j+1, jb, j-1, -one, a( j, 1 ), lda, a( 1, j ), lda, one, a( j, j ), lda )

             CALL dgetf2( m-j+1, jb, a( j, j ), lda, ipiv( j ), iinfo )

             IF( info.EQ.0 .AND. iinfo.GT.0 ) info = iinfo + j - 1
             DO 10 i = j, min( m, j+jb-1 )
                ipiv( i ) = j - 1 + ipiv( i )
    10       CONTINUE

             CALL dlaswp( j-1, a, lda, j, j+jb-1, ipiv, 1 )

             IF ( j+jb.LE.n ) THEN

                CALL dlaswp( n-j-jb+1, a( 1, j+jb ), lda, j, j+jb-1, ipiv, 1 )

                CALL dgemm( 'No transpose', 'No transpose', jb, n-j-jb+1, j-1, -one, a( j, 1 ), lda, a( 1, j+jb ), lda, one, a( j, j+jb ), lda )

                CALL dtrsm( 'Left', 'Lower', 'No transpose', 'Unit', jb, n-j-jb+1, one, a( j, j ), lda, a( j, j+jb ), lda )
             END IF

    20    CONTINUE

       END IF
       RETURN

       END SUBROUTINE dgetrf

       SUBROUTINE dtrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

       DOUBLE PRECISION ALPHA
       INTEGER LDA,LDB,M,N
       CHARACTER DIAG,SIDE,TRANSA,UPLO

       DOUBLE PRECISION A(LDA,*),B(LDB,*)

       !LOGICAL LSAME

       INTRINSIC max

       DOUBLE PRECISION TEMP
       INTEGER I,INFO,J,K,NROWA
       LOGICAL LSIDE,NOUNIT,UPPER

       DOUBLE PRECISION ONE,ZERO
       parameter(one=1.0d+0,zero=0.0d+0)

       lside = lsame(side,'L')
       IF (lside) THEN
           nrowa = m
       ELSE
           nrowa = n
       END IF
       nounit = lsame(diag,'N')
       upper = lsame(uplo,'U')

       info = 0
       IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
           info = 1
       ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
           info = 2
       ELSE IF ((.NOT.lsame(transa,'N')) .AND. + (.NOT.lsame(transa,'T')) .AND. + (.NOT.lsame(transa,'C'))) THEN
           info = 3
       ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
           info = 4
       ELSE IF (m.LT.0) THEN
           info = 5
       ELSE IF (n.LT.0) THEN
           info = 6
       ELSE IF (lda.LT.max(1,nrowa)) THEN
           info = 9
       ELSE IF (ldb.LT.max(1,m)) THEN
           info = 11
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('DTRSM ',info)
           RETURN
       END IF

       IF (m.EQ.0 .OR. n.EQ.0) RETURN

       IF (alpha.EQ.zero) THEN
           DO 20 j = 1,n
               DO 10 i = 1,m
                   b(i,j) = zero
    10         CONTINUE
    20     CONTINUE
           RETURN
       END IF

       IF (lside) THEN
           IF (lsame(transa,'N')) THEN

               IF (upper) THEN
                   DO 60 j = 1,n
                       IF (alpha.NE.one) THEN
                           DO 30 i = 1,m
                               b(i,j) = alpha*b(i,j)
    30                     CONTINUE
                       END IF
                       DO 50 k = m,1,-1
                           IF (b(k,j).NE.zero) THEN
                               IF (nounit) b(k,j) = b(k,j)/a(k,k)
                               DO 40 i = 1,k - 1
                                   b(i,j) = b(i,j) - b(k,j)*a(i,k)
    40                         CONTINUE
                           END IF
    50                 CONTINUE
    60             CONTINUE
               ELSE
                   DO 100 j = 1,n
                       IF (alpha.NE.one) THEN
                           DO 70 i = 1,m
                               b(i,j) = alpha*b(i,j)
    70                     CONTINUE
                       END IF
                       DO 90 k = 1,m
                           IF (b(k,j).NE.zero) THEN
                               IF (nounit) b(k,j) = b(k,j)/a(k,k)
                               DO 80 i = k + 1,m
                                   b(i,j) = b(i,j) - b(k,j)*a(i,k)
    80                         CONTINUE
                           END IF
    90                 CONTINUE
   100             CONTINUE
               END IF
           ELSE

               IF (upper) THEN
                   DO 130 j = 1,n
                       DO 120 i = 1,m
                           temp = alpha*b(i,j)
                           DO 110 k = 1,i - 1
                               temp = temp - a(k,i)*b(k,j)
   110                     CONTINUE
                           IF (nounit) temp = temp/a(i,i)
                           b(i,j) = temp
   120                 CONTINUE
   130             CONTINUE
               ELSE
                   DO 160 j = 1,n
                       DO 150 i = m,1,-1
                           temp = alpha*b(i,j)
                           DO 140 k = i + 1,m
                               temp = temp - a(k,i)*b(k,j)
   140                     CONTINUE
                           IF (nounit) temp = temp/a(i,i)
                           b(i,j) = temp
   150                 CONTINUE
   160             CONTINUE
               END IF
           END IF
       ELSE
           IF (lsame(transa,'N')) THEN

               IF (upper) THEN
                   DO 210 j = 1,n
                       IF (alpha.NE.one) THEN
                           DO 170 i = 1,m
                               b(i,j) = alpha*b(i,j)
   170                     CONTINUE
                       END IF
                       DO 190 k = 1,j - 1
                           IF (a(k,j).NE.zero) THEN
                               DO 180 i = 1,m
                                   b(i,j) = b(i,j) - a(k,j)*b(i,k)
   180                         CONTINUE
                           END IF
   190                 CONTINUE
                       IF (nounit) THEN
                           temp = one/a(j,j)
                           DO 200 i = 1,m
                               b(i,j) = temp*b(i,j)
   200                     CONTINUE
                       END IF
   210             CONTINUE
               ELSE
                   DO 260 j = n,1,-1
                       IF (alpha.NE.one) THEN
                           DO 220 i = 1,m
                               b(i,j) = alpha*b(i,j)
   220                     CONTINUE
                       END IF
                       DO 240 k = j + 1,n
                           IF (a(k,j).NE.zero) THEN
                               DO 230 i = 1,m
                                   b(i,j) = b(i,j) - a(k,j)*b(i,k)
   230                         CONTINUE
                           END IF
   240                 CONTINUE
                       IF (nounit) THEN
                           temp = one/a(j,j)
                           DO 250 i = 1,m
                               b(i,j) = temp*b(i,j)
   250                     CONTINUE
                       END IF
   260             CONTINUE
               END IF
           ELSE

               IF (upper) THEN
                   DO 310 k = n,1,-1
                       IF (nounit) THEN
                           temp = one/a(k,k)
                           DO 270 i = 1,m
                               b(i,k) = temp*b(i,k)
   270                     CONTINUE
                       END IF
                       DO 290 j = 1,k - 1
                           IF (a(j,k).NE.zero) THEN
                               temp = a(j,k)
                               DO 280 i = 1,m
                                   b(i,j) = b(i,j) - temp*b(i,k)
   280                         CONTINUE
                           END IF
   290                 CONTINUE
                       IF (alpha.NE.one) THEN
                           DO 300 i = 1,m
                               b(i,k) = alpha*b(i,k)
   300                     CONTINUE
                       END IF
   310             CONTINUE
               ELSE
                   DO 360 k = 1,n
                       IF (nounit) THEN
                           temp = one/a(k,k)
                           DO 320 i = 1,m
                               b(i,k) = temp*b(i,k)
   320                     CONTINUE
                       END IF
                       DO 340 j = k + 1,n
                           IF (a(j,k).NE.zero) THEN
                               temp = a(j,k)
                               DO 330 i = 1,m
                                   b(i,j) = b(i,j) - temp*b(i,k)
   330                         CONTINUE
                           END IF
   340                 CONTINUE
                       IF (alpha.NE.one) THEN
                           DO 350 i = 1,m
                               b(i,k) = alpha*b(i,k)
   350                     CONTINUE
                       END IF
   360             CONTINUE
               END IF
           END IF
       END IF

       RETURN

       END SUBROUTINE dtrsm

       SUBROUTINE dgetf2( M, N, A, LDA, IPIV, INFO )

       INTEGER            INFO, LDA, M, N

       INTEGER            IPIV( * )
       DOUBLE PRECISION   A( LDA, * )

       DOUBLE PRECISION   ONE, ZERO
       parameter( one = 1.0d+0, zero = 0.0d+0 )

       DOUBLE PRECISION   SFMIN
       INTEGER            I, J, JP

       !DOUBLE PRECISION   DLAMCH
       !INTEGER            IDAMAX

       INTRINSIC          max, min

       info = 0
       IF( m.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = -4
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'DGETF2', -info )
          RETURN
       END IF

       IF( m.EQ.0 .OR. n.EQ.0 ) RETURN

       sfmin = dlamch('S')

       DO 10 j = 1, min( m, n )

          jp = j - 1 + idamax( m-j+1, a( j, j ), 1 )
          ipiv( j ) = jp
          IF( a( jp, j ).NE.zero ) THEN

             IF( jp.NE.j ) CALL dswap( n, a( j, 1 ), lda, a( jp, 1 ), lda )

             IF( j.LT.m ) THEN
                IF( abs(a( j, j )) .GE. sfmin ) THEN
                   CALL dscal( m-j, one / a( j, j ), a( j+1, j ), 1 )
                ELSE
                  DO 20 i = 1, m-j
                     a( j+i, j ) = a( j+i, j ) / a( j, j )
    20            CONTINUE
                END IF
             END IF

          ELSE IF( info.EQ.0 ) THEN

             info = j
          END IF

          IF( j.LT.min( m, n ) ) THEN

             CALL dger( m-j, n-j, -one, a( j+1, j ), 1, a( j, j+1 ), lda, a( j+1, j+1 ), lda )
          END IF
    10 CONTINUE
       RETURN

       END SUBROUTINE dgetf2

       SUBROUTINE dlaswp( N, A, LDA, K1, K2, IPIV, INCX )

       INTEGER            INCX, K1, K2, LDA, N

       INTEGER            IPIV( * )
       DOUBLE PRECISION   A( LDA, * )

       INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
       DOUBLE PRECISION   TEMP

       IF( incx.GT.0 ) THEN
          ix0 = k1
          i1 = k1
          i2 = k2
          inc = 1
       ELSE IF( incx.LT.0 ) THEN
          ix0 = k1 + ( k1-k2 )*incx
          i1 = k2
          i2 = k1
          inc = -1
       ELSE
          RETURN
       END IF

       n32 = ( n / 32 )*32
       IF( n32.NE.0 ) THEN
          DO 30 j = 1, n32, 32
             ix = ix0
             DO 20 i = i1, i2, inc
                ip = ipiv( ix )
                IF( ip.NE.i ) THEN
                   DO 10 k = j, j + 31
                      temp = a( i, k )
                      a( i, k ) = a( ip, k )
                      a( ip, k ) = temp
    10             CONTINUE
                END IF
                ix = ix + incx
    20       CONTINUE
    30    CONTINUE
       END IF
       IF( n32.NE.n ) THEN
          n32 = n32 + 1
          ix = ix0
          DO 50 i = i1, i2, inc
             ip = ipiv( ix )
             IF( ip.NE.i ) THEN
                DO 40 k = n32, n
                   temp = a( i, k )
                   a( i, k ) = a( ip, k )
                   a( ip, k ) = temp
    40          CONTINUE
             END IF
             ix = ix + incx
    50    CONTINUE
       END IF

       RETURN

       END SUBROUTINE dlaswp

       SUBROUTINE dger(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)

       DOUBLE PRECISION ALPHA
       INTEGER INCX,INCY,LDA,M,N

       DOUBLE PRECISION A(LDA,*),X(*),Y(*)

       DOUBLE PRECISION ZERO
       parameter(zero=0.0d+0)

       DOUBLE PRECISION TEMP
       INTEGER I,INFO,IX,J,JY,KX

       INTRINSIC max

       info = 0
       IF (m.LT.0) THEN
           info = 1
       ELSE IF (n.LT.0) THEN
           info = 2
       ELSE IF (incx.EQ.0) THEN
           info = 5
       ELSE IF (incy.EQ.0) THEN
           info = 7
       ELSE IF (lda.LT.max(1,m)) THEN
           info = 9
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('DGER  ',info)
           RETURN
       END IF

       IF ((m.EQ.0) .OR. (n.EQ.0) .OR. (alpha.EQ.zero)) RETURN

       IF (incy.GT.0) THEN
           jy = 1
       ELSE
           jy = 1 - (n-1)*incy
       END IF
       IF (incx.EQ.1) THEN
           DO 20 j = 1,n
               IF (y(jy).NE.zero) THEN
                   temp = alpha*y(jy)
                   DO 10 i = 1,m
                       a(i,j) = a(i,j) + x(i)*temp
    10             CONTINUE
               END IF
               jy = jy + incy
    20     CONTINUE
       ELSE
           IF (incx.GT.0) THEN
               kx = 1
           ELSE
               kx = 1 - (m-1)*incx
           END IF
           DO 40 j = 1,n
               IF (y(jy).NE.zero) THEN
                   temp = alpha*y(jy)
                   ix = kx
                   DO 30 i = 1,m
                       a(i,j) = a(i,j) + x(ix)*temp
                       ix = ix + incx
    30             CONTINUE
               END IF
               jy = jy + incy
    40     CONTINUE
       END IF

       RETURN

       END SUBROUTINE dger

       SUBROUTINE dscal(N,DA,DX,INCX)

       DOUBLE PRECISION DA
       INTEGER INCX,N

       DOUBLE PRECISION DX(*)

       INTEGER I,M,MP1,NINCX

       INTRINSIC mod

       IF (n.LE.0 .OR. incx.LE.0) RETURN
       IF (incx.EQ.1) THEN

          m = mod(n,5)
          IF (m.NE.0) THEN
             DO i = 1,m
                dx(i) = da*dx(i)
             END DO
             IF (n.LT.5) RETURN
          END IF
          mp1 = m + 1
          DO i = mp1,n,5
             dx(i) = da*dx(i)
             dx(i+1) = da*dx(i+1)
             dx(i+2) = da*dx(i+2)
             dx(i+3) = da*dx(i+3)
             dx(i+4) = da*dx(i+4)
          END DO
       ELSE

          nincx = n*incx
          DO i = 1,nincx,incx
             dx(i) = da*dx(i)
          END DO
       END IF
       RETURN

       END SUBROUTINE dscal

       SUBROUTINE dswap(N,DX,INCX,DY,INCY)

       INTEGER INCX,INCY,N

       DOUBLE PRECISION DX(*),DY(*)

       DOUBLE PRECISION DTEMP
       INTEGER I,IX,IY,M,MP1

       INTRINSIC mod

       IF (n.LE.0) RETURN
       IF (incx.EQ.1 .AND. incy.EQ.1) THEN

          m = mod(n,3)
          IF (m.NE.0) THEN
             DO i = 1,m
                dtemp = dx(i)
                dx(i) = dy(i)
                dy(i) = dtemp
             END DO
             IF (n.LT.3) RETURN
          END IF
          mp1 = m + 1
          DO i = mp1,n,3
             dtemp = dx(i)
             dx(i) = dy(i)
             dy(i) = dtemp
             dtemp = dx(i+1)
             dx(i+1) = dy(i+1)
             dy(i+1) = dtemp
             dtemp = dx(i+2)
             dx(i+2) = dy(i+2)
             dy(i+2) = dtemp
          END DO
       ELSE

          ix = 1
          iy = 1
          IF (incx.LT.0) ix = (-n+1)*incx + 1
          IF (incy.LT.0) iy = (-n+1)*incy + 1
          DO i = 1,n
             dtemp = dx(ix)
             dx(ix) = dy(iy)
             dy(iy) = dtemp
             ix = ix + incx
             iy = iy + incy
          END DO
       END IF
       RETURN

       END SUBROUTINE dswap

       DOUBLE PRECISION FUNCTION dlamch( CMACH )

       CHARACTER          cmach

       DOUBLE PRECISION   one, zero
       parameter( one = 1.0d+0, zero = 0.0d+0 )

       DOUBLE PRECISION   rnd, eps, sfmin, small, rmach

       !LOGICAL            lsame

       INTRINSIC          digits, epsilon, huge, maxexponent, minexponent, radix, tiny

       rnd = one

       IF( one.EQ.rnd ) THEN
          eps = epsilon(zero) * 0.5
       ELSE
          eps = epsilon(zero)
       END IF

       IF( lsame( cmach, 'E' ) ) THEN
          rmach = eps
       ELSE IF( lsame( cmach, 'S' ) ) THEN
          sfmin = tiny(zero)
          small = one / huge(zero)
          IF( small.GE.sfmin ) THEN
             sfmin = small*( one+eps )
          END IF
          rmach = sfmin
       ELSE IF( lsame( cmach, 'B' ) ) THEN
          rmach = radix(zero)
       ELSE IF( lsame( cmach, 'P' ) ) THEN
          rmach = eps * radix(zero)
       ELSE IF( lsame( cmach, 'N' ) ) THEN
          rmach = digits(zero)
       ELSE IF( lsame( cmach, 'R' ) ) THEN
          rmach = rnd
       ELSE IF( lsame( cmach, 'M' ) ) THEN
          rmach = minexponent(zero)
       ELSE IF( lsame( cmach, 'U' ) ) THEN
          rmach = tiny(zero)
       ELSE IF( lsame( cmach, 'L' ) ) THEN
          rmach = maxexponent(zero)
       ELSE IF( lsame( cmach, 'O' ) ) THEN
          rmach = huge(zero)
       ELSE
          rmach = zero
       END IF

       dlamch = rmach
       RETURN

       END FUNCTION dlamch

       DOUBLE PRECISION FUNCTION dlamc3( A, B )

       DOUBLE PRECISION   a, b

       dlamc3 = a + b

       RETURN

       END FUNCTION dlamc3

       INTEGER FUNCTION ilaenv( ISPEC, NAME, OPTS, N1, N2, N3, N4 )

       CHARACTER*( * )    name, opts
       INTEGER            ispec, n1, n2, n3, n4

       INTEGER            i, ic, iz, nb, nbmin, nx
       LOGICAL            cname, sname, twostage
       CHARACTER          c1*1, c2*2, c4*2, c3*3, subnam*16

       INTRINSIC          char, ichar, int, min, real

       !INTEGER            ieeeck, iparmq, iparam2stage

       GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, 130, 140, 150, 160, 160, 160, 160, 160, 160)ispec

       ilaenv = -1
       RETURN

    10 CONTINUE

       ilaenv = 1
       subnam = name
       ic = ichar( subnam( 1: 1 ) )
       iz = ichar( 'Z' )
       IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN

          IF( ic.GE.97 .AND. ic.LE.122 ) THEN
             subnam( 1: 1 ) = char( ic-32 )
             DO 20 i = 2, 6
                ic = ichar( subnam( i: i ) )
                IF( ic.GE.97 .AND. ic.LE.122 ) subnam( i: i ) = char( ic-32 )
    20       CONTINUE
          END IF

       ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN

          IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. ( ic.GE.145 .AND. ic.LE.153 ) .OR. ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
             subnam( 1: 1 ) = char( ic+64 )
             DO 30 i = 2, 6
                ic = ichar( subnam( i: i ) )
                IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. ( ic.GE.145 .AND. ic.LE.153 ) .OR. ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i:i ) = char( ic+64 )
    30       CONTINUE
          END IF

       ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN

          IF( ic.GE.225 .AND. ic.LE.250 ) THEN
             subnam( 1: 1 ) = char( ic-32 )
             DO 40 i = 2, 6
                ic = ichar( subnam( i: i ) )
                IF( ic.GE.225 .AND. ic.LE.250 ) subnam( i: i ) = char( ic-32 )
    40       CONTINUE
          END IF
       END IF

       c1 = subnam( 1: 1 )
       sname = c1.EQ.'S' .OR. c1.EQ.'D'
       cname = c1.EQ.'C' .OR. c1.EQ.'Z'
       IF( .NOT.( cname .OR. sname ) ) RETURN
       c2 = subnam( 2: 3 )
       c3 = subnam( 4: 6 )
       c4 = c3( 2: 3 )
       twostage = len( subnam ).GE.11 .AND. subnam( 11: 11 ).EQ.'2'

       GO TO ( 50, 60, 70 )ispec

    50 CONTINUE

       nb = 1

       IF( subnam(2:6).EQ.'LAORH' ) THEN

          IF( sname ) THEN
              nb = 32
          ELSE
              nb = 32
          END IF
       ELSE IF( c2.EQ.'GE' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( sname ) THEN
                nb = 64
             ELSE
                nb = 64
             END IF
          ELSE IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ.'QLF' ) THEN
             IF( sname ) THEN
                nb = 32
             ELSE
                nb = 32
             END IF
          ELSE IF( c3.EQ.'QR ') THEN
             IF( n3 .EQ. 1) THEN
                IF( sname ) THEN
                   IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                      nb = n1
                   ELSE
                      nb = 32768/n2
                   END IF
                ELSE
                   IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                      nb = n1
                   ELSE
                      nb = 32768/n2
                   END IF
                END IF
             ELSE
                IF( sname ) THEN
                   nb = 1
                ELSE
                   nb = 1
                END IF
             END IF
          ELSE IF( c3.EQ.'LQ ') THEN
             IF( n3 .EQ. 2) THEN
                IF( sname ) THEN
                   IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                      nb = n1
                   ELSE
                      nb = 32768/n2
                   END IF
                ELSE
                   IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                      nb = n1
                   ELSE
                      nb = 32768/n2
                   END IF
                END IF
             ELSE
                IF( sname ) THEN
                   nb = 1
                ELSE
                   nb = 1
                END IF
             END IF
          ELSE IF( c3.EQ.'HRD' ) THEN
             IF( sname ) THEN
                nb = 32
             ELSE
                nb = 32
             END IF
          ELSE IF( c3.EQ.'BRD' ) THEN
             IF( sname ) THEN
                nb = 32
             ELSE
                nb = 32
             END IF
          ELSE IF( c3.EQ.'TRI' ) THEN
             IF( sname ) THEN
                nb = 64
             ELSE
                nb = 64
             END IF
          END IF
       ELSE IF( c2.EQ.'PO' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( sname ) THEN
                nb = 64
             ELSE
                nb = 64
             END IF
          END IF
       ELSE IF( c2.EQ.'SY' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( sname ) THEN
                IF( twostage ) THEN
                   nb = 192
                ELSE
                   nb = 64
                END IF
             ELSE
                IF( twostage ) THEN
                   nb = 192
                ELSE
                   nb = 64
                END IF
             END IF
          ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
             nb = 32
          ELSE IF( sname .AND. c3.EQ.'GST' ) THEN
             nb = 64
          END IF
       ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( twostage ) THEN
                nb = 192
             ELSE
                nb = 64
             END IF
          ELSE IF( c3.EQ.'TRD' ) THEN
             nb = 32
          ELSE IF( c3.EQ.'GST' ) THEN
             nb = 64
          END IF
       ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
          IF( c3( 1: 1 ).EQ.'G' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nb = 32
             END IF
          ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )THEN
                nb = 32
             END IF
          END IF
       ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
          IF( c3( 1: 1 ).EQ.'G' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nb = 32
             END IF
          ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nb = 32
             END IF
          END IF
       ELSE IF( c2.EQ.'GB' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( sname ) THEN
                IF( n4.LE.64 ) THEN
                   nb = 1
                ELSE
                   nb = 32
                END IF
             ELSE
                IF( n4.LE.64 ) THEN
                   nb = 1
                ELSE
                   nb = 32
                END IF
             END IF
          END IF
       ELSE IF( c2.EQ.'PB' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( sname ) THEN
                IF( n2.LE.64 ) THEN
                   nb = 1
                ELSE
                   nb = 32
                END IF
             ELSE
                IF( n2.LE.64 ) THEN
                   nb = 1
                ELSE
                   nb = 32
                END IF
             END IF
          END IF
       ELSE IF( c2.EQ.'TR' ) THEN
          IF( c3.EQ.'TRI' ) THEN
             IF( sname ) THEN
                nb = 64
             ELSE
                nb = 64
             END IF
          ELSE IF ( c3.EQ.'EVC' ) THEN
             IF( sname ) THEN
                nb = 64
             ELSE
                nb = 64
             END IF
          END IF
       ELSE IF( c2.EQ.'LA' ) THEN
          IF( c3.EQ.'UUM' ) THEN
             IF( sname ) THEN
                nb = 64
             ELSE
                nb = 64
             END IF
          END IF
       ELSE IF( sname .AND. c2.EQ.'ST' ) THEN
          IF( c3.EQ.'EBZ' ) THEN
             nb = 1
          END IF
       ELSE IF( c2.EQ.'GG' ) THEN
          nb = 32
          IF( c3.EQ.'HD3' ) THEN
             IF( sname ) THEN
                nb = 32
             ELSE
                nb = 32
             END IF
          END IF
       END IF
       ilaenv = nb
       RETURN

    60 CONTINUE

       nbmin = 2
       IF( c2.EQ.'GE' ) THEN
          IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ.'QLF' ) THEN
             IF( sname ) THEN
                nbmin = 2
             ELSE
                nbmin = 2
             END IF
          ELSE IF( c3.EQ.'HRD' ) THEN
             IF( sname ) THEN
                nbmin = 2
             ELSE
                nbmin = 2
             END IF
          ELSE IF( c3.EQ.'BRD' ) THEN
             IF( sname ) THEN
                nbmin = 2
             ELSE
                nbmin = 2
             END IF
          ELSE IF( c3.EQ.'TRI' ) THEN
             IF( sname ) THEN
                nbmin = 2
             ELSE
                nbmin = 2
             END IF
          END IF
       ELSE IF( c2.EQ.'SY' ) THEN
          IF( c3.EQ.'TRF' ) THEN
             IF( sname ) THEN
                nbmin = 8
             ELSE
                nbmin = 8
             END IF
          ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
             nbmin = 2
          END IF
       ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
          IF( c3.EQ.'TRD' ) THEN
             nbmin = 2
          END IF
       ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
          IF( c3( 1: 1 ).EQ.'G' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nbmin = 2
             END IF
          ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nbmin = 2
             END IF
          END IF
       ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
          IF( c3( 1: 1 ).EQ.'G' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nbmin = 2
             END IF
          ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nbmin = 2
             END IF
          END IF
       ELSE IF( c2.EQ.'GG' ) THEN
          nbmin = 2
          IF( c3.EQ.'HD3' ) THEN
             nbmin = 2
          END IF
       END IF
       ilaenv = nbmin
       RETURN

    70 CONTINUE

       nx = 0
       IF( c2.EQ.'GE' ) THEN
          IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ.'QLF' ) THEN
             IF( sname ) THEN
                nx = 128
             ELSE
                nx = 128
             END IF
          ELSE IF( c3.EQ.'HRD' ) THEN
             IF( sname ) THEN
                nx = 128
             ELSE
                nx = 128
             END IF
          ELSE IF( c3.EQ.'BRD' ) THEN
             IF( sname ) THEN
                nx = 128
             ELSE
                nx = 128
             END IF
          END IF
       ELSE IF( c2.EQ.'SY' ) THEN
          IF( sname .AND. c3.EQ.'TRD' ) THEN
             nx = 32
          END IF
       ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
          IF( c3.EQ.'TRD' ) THEN
             nx = 32
          END IF
       ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
          IF( c3( 1: 1 ).EQ.'G' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) THEN
                nx = 128
             END IF
          END IF
       ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
          IF( c3( 1: 1 ).EQ.'G' ) THEN
             IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )THEN
                nx = 128
             END IF
          END IF
       ELSE IF( c2.EQ.'GG' ) THEN
          nx = 128
          IF( c3.EQ.'HD3' ) THEN
             nx = 128
          END IF
       END IF
       ilaenv = nx
       RETURN

    80 CONTINUE

       ilaenv = 6
       RETURN

    90 CONTINUE

       ilaenv = 2
       RETURN

   100 CONTINUE

       ilaenv = int( real( min( n1, n2 ) )*1.6e0 )
       RETURN

   110 CONTINUE

       ilaenv = 1
       RETURN

   120 CONTINUE

       ilaenv = 50
       RETURN

   130 CONTINUE

       ilaenv = 25
       RETURN

   140 CONTINUE

       ilaenv = 1
       IF( ilaenv.EQ.1 ) THEN
          ilaenv = ieeeck( 1, 0.0, 1.0 )
       END IF
       RETURN

   150 CONTINUE

       ilaenv = 1
       IF( ilaenv.EQ.1 ) THEN
          ilaenv = ieeeck( 0, 0.0, 1.0 )
       END IF
       RETURN

   160 CONTINUE

       ilaenv = iparmq( ispec, name, opts, n1, n2, n3, n4 )
       RETURN

       END FUNCTION ilaenv

       INTEGER FUNCTION idamax(N,DX,INCX)

       INTEGER incx,n

       DOUBLE PRECISION dx(*)

       DOUBLE PRECISION dmax
       INTEGER i,ix

       INTRINSIC dabs

       idamax = 0
       IF (n.LT.1 .OR. incx.LE.0) RETURN
       idamax = 1
       IF (n.EQ.1) RETURN
       IF (incx.EQ.1) THEN

          dmax = dabs(dx(1))
          DO i = 2,n
             IF (dabs(dx(i)).GT.dmax) THEN
                idamax = i
                dmax = dabs(dx(i))
             END IF
          END DO
       ELSE

          ix = 1
          dmax = dabs(dx(1))
          ix = ix + incx
          DO i = 2,n
             IF (dabs(dx(ix)).GT.dmax) THEN
                idamax = i
                dmax = dabs(dx(ix))
             END IF
             ix = ix + incx
          END DO
       END IF
       RETURN

       END FUNCTION idamax

       INTEGER FUNCTION ieeeck( ISPEC, ZERO, ONE )

       INTEGER            ispec
       REAL               one, zero

       REAL               nan1, nan2, nan3, nan4, nan5, nan6, neginf, negzro, newzro, posinf

       ieeeck = 1

       posinf = one / zero
       IF( posinf.LE.one ) THEN
          ieeeck = 0
          RETURN
       END IF

       neginf = -one / zero
       IF( neginf.GE.zero ) THEN
          ieeeck = 0
          RETURN
       END IF

       negzro = one / ( neginf+one )
       IF( negzro.NE.zero ) THEN
          ieeeck = 0
          RETURN
       END IF

       neginf = one / negzro
       IF( neginf.GE.zero ) THEN
          ieeeck = 0
          RETURN
       END IF

       newzro = negzro + zero
       IF( newzro.NE.zero ) THEN
          ieeeck = 0
          RETURN
       END IF

       posinf = one / newzro
       IF( posinf.LE.one ) THEN
          ieeeck = 0
          RETURN
       END IF

       neginf = neginf*posinf
       IF( neginf.GE.zero ) THEN
          ieeeck = 0
          RETURN
       END IF

       posinf = posinf*posinf
       IF( posinf.LE.one ) THEN
          ieeeck = 0
          RETURN
       END IF

       IF( ispec.EQ.0 ) RETURN

       nan1 = posinf + neginf

       nan2 = posinf / neginf

       nan3 = posinf / posinf

       nan4 = posinf*zero

       nan5 = neginf*negzro

       nan6 = nan5*zero

       IF( nan1.EQ.nan1 ) THEN
          ieeeck = 0
          RETURN
       END IF

       IF( nan2.EQ.nan2 ) THEN
          ieeeck = 0
          RETURN
       END IF

       IF( nan3.EQ.nan3 ) THEN
          ieeeck = 0
          RETURN
       END IF

       IF( nan4.EQ.nan4 ) THEN
          ieeeck = 0
          RETURN
       END IF

       IF( nan5.EQ.nan5 ) THEN
          ieeeck = 0
          RETURN
       END IF

       IF( nan6.EQ.nan6 ) THEN
          ieeeck = 0
          RETURN
       END IF

       RETURN
       END function ieeeck

       INTEGER FUNCTION iparmq( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )

       INTEGER            ihi, ilo, ispec, lwork, n
       CHARACTER          name*( * ), opts*( * )

       INTEGER            inmin, inwin, inibl, ishfts, iacc22, icost
       parameter( inmin = 12, inwin = 13, inibl = 14, ishfts = 15, iacc22 = 16, icost = 17 )
       INTEGER            nmin, k22min, kacmin, nibble, knwswp, rcost
       parameter( nmin = 75, k22min = 14, kacmin = 14, nibble = 14, knwswp = 500, rcost = 10 )
       REAL               two
       parameter( two = 2.0 )

       INTEGER            nh, ns
       INTEGER            i, ic, iz
       CHARACTER          subnam*6

       INTRINSIC          log, max, mod, nint, real

       IF( ( ispec.EQ.ishfts ) .OR. ( ispec.EQ.inwin ) .OR. ( ispec.EQ.iacc22 ) ) THEN

          nh = ihi - ilo + 1
          ns = 2
          IF( nh.GE.30 ) ns = 4
          IF( nh.GE.60 ) ns = 10
          IF( nh.GE.150 ) ns = max( 10, nh / nint( log( real( nh ) ) / log( two ) ) )
          IF( nh.GE.590 ) ns = 64
          IF( nh.GE.3000 ) ns = 128
          IF( nh.GE.6000 ) ns = 256
          ns = max( 2, ns-mod( ns, 2 ) )
       END IF

       IF( ispec.EQ.inmin ) THEN

          iparmq = nmin

       ELSE IF( ispec.EQ.inibl ) THEN

          iparmq = nibble

       ELSE IF( ispec.EQ.ishfts ) THEN

          iparmq = ns

       ELSE IF( ispec.EQ.inwin ) THEN

          IF( nh.LE.knwswp ) THEN
             iparmq = ns
          ELSE
             iparmq = 3*ns / 2
          END IF

       ELSE IF( ispec.EQ.iacc22 ) THEN

          iparmq = 0
          subnam = name
          ic = ichar( subnam( 1: 1 ) )
          iz = ichar( 'Z' )
          IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN

             IF( ic.GE.97 .AND. ic.LE.122 ) THEN
                subnam( 1: 1 ) = char( ic-32 )
                DO i = 2, 6
                   ic = ichar( subnam( i: i ) )
                   IF( ic.GE.97 .AND. ic.LE.122 ) subnam( i: i ) = char( ic-32 )
                END DO
             END IF

          ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN

             IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. ( ic.GE.145 .AND. ic.LE.153 ) .OR. ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
                subnam( 1: 1 ) = char( ic+64 )
                DO i = 2, 6
                   ic = ichar( subnam( i: i ) )
                   IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. ( ic.GE.145 .AND. ic.LE.153 ) .OR. ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i:i ) = char( ic+64 )
                END DO
             END IF

          ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN

             IF( ic.GE.225 .AND. ic.LE.250 ) THEN
                subnam( 1: 1 ) = char( ic-32 )
                DO i = 2, 6
                   ic = ichar( subnam( i: i ) )
                   IF( ic.GE.225 .AND. ic.LE.250 ) subnam( i: i ) = char( ic-32 )
                END DO
             END IF
          END IF

          IF( subnam( 2:6 ).EQ.'GGHRD' .OR. subnam( 2:6 ).EQ.'GGHD3' ) THEN
             iparmq = 1
             IF( nh.GE.k22min ) iparmq = 2
          ELSE IF ( subnam( 4:6 ).EQ.'EXC' ) THEN
             IF( nh.GE.kacmin ) iparmq = 1
             IF( nh.GE.k22min ) iparmq = 2
          ELSE IF ( subnam( 2:6 ).EQ.'HSEQR' .OR. subnam( 2:5 ).EQ.'LAQR' ) THEN
             IF( ns.GE.kacmin ) iparmq = 1
             IF( ns.GE.k22min ) iparmq = 2
          END IF

       ELSE IF( ispec.EQ.icost ) THEN

          iparmq = rcost
       ELSE
          iparmq = -1
       END IF

       END FUNCTION iparmq

       INTEGER FUNCTION iparam2stage( ISPEC, NAME, OPTS, NI, NBI, IBI, NXI )

       IMPLICIT NONE

       CHARACTER*( * )    name, opts
       INTEGER            ispec, ni, nbi, ibi, nxi

       INTEGER            i, ic, iz, kd, ib, lhous, lwork, nthreads, factoptnb, qroptnb, lqoptnb
       LOGICAL            rprec, cprec
       CHARACTER          prec*1, algo*3, stag*5, subnam*12, vect*1

       INTRINSIC          char, ichar, max

       !INTEGER            ilaenv

       IF( (ispec.LT.17).OR.(ispec.GT.21) ) THEN
           iparam2stage = -1
           RETURN
       ENDIF

       nthreads = 1

       IF( ispec .NE. 19 ) THEN
          iparam2stage = -1
          subnam = name
          ic = ichar( subnam( 1: 1 ) )
          iz = ichar( 'Z' )
          IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
             IF( ic.GE.97 .AND. ic.LE.122 ) THEN
                subnam( 1: 1 ) = char( ic-32 )
                DO 100 i = 2, 12
                   ic = ichar( subnam( i: i ) )
                   IF( ic.GE.97 .AND. ic.LE.122 )subnam( i: i ) = char( ic-32 )
   100          CONTINUE
             END IF
          ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
             IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. ( ic.GE.145 .AND. ic.LE.153 ) .OR. ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
                subnam( 1: 1 ) = char( ic+64 )
                DO 110 i = 2, 12
                   ic = ichar( subnam( i: i ) )
                   IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. ( ic.GE.145 .AND. ic.LE.153 ) .OR. ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i:i ) = char( ic+64 )
   110          CONTINUE
             END IF

          ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN

             IF( ic.GE.225 .AND. ic.LE.250 ) THEN
                subnam( 1: 1 ) = char( ic-32 )
                DO 120 i = 2, 12
                  ic = ichar( subnam( i: i ) )
                  IF( ic.GE.225 .AND. ic.LE.250 ) subnam( i: i ) = char( ic-32 )
   120          CONTINUE
             END IF
          END IF

          prec  = subnam( 1: 1 )
          algo  = subnam( 4: 6 )
          stag  = subnam( 8:12 )
          rprec = prec.EQ.'S' .OR. prec.EQ.'D'
          cprec = prec.EQ.'C' .OR. prec.EQ.'Z'

          IF( .NOT.( rprec .OR. cprec ) ) THEN
              iparam2stage = -1
              RETURN
          ENDIF
       ENDIF

       IF (( ispec .EQ. 17 ) .OR. ( ispec .EQ. 18 )) THEN

          IF( nthreads.GT.4 ) THEN
             IF( cprec ) THEN
                kd = 128
                ib = 32
             ELSE
                kd = 160
                ib = 40
             ENDIF
          ELSE IF( nthreads.GT.1 ) THEN
             IF( cprec ) THEN
                kd = 64
                ib = 32
             ELSE
                kd = 64
                ib = 32
             ENDIF
          ELSE
             IF( cprec ) THEN
                kd = 16
                ib = 16
             ELSE
                kd = 32
                ib = 16
             ENDIF
          ENDIF
          IF( ispec.EQ.17 ) iparam2stage = kd
          IF( ispec.EQ.18 ) iparam2stage = ib

       ELSE IF ( ispec .EQ. 19 ) THEN

          vect  = opts(1:1)
          IF( vect.EQ.'N' ) THEN
             lhous = max( 1, 4*ni )
          ELSE
             lhous = max( 1, 4*ni ) + ibi
          ENDIF
          IF( lhous.GE.0 ) THEN
             iparam2stage = lhous
          ELSE
             iparam2stage = -1
          ENDIF

       ELSE IF ( ispec .EQ. 20 ) THEN
          lwork        = -1
          subnam(1:1)  = prec
          subnam(2:6)  = 'GEQRF'
          qroptnb      = ilaenv( 1, subnam, ' ', ni, nbi, -1, -1 )
          subnam(2:6)  = 'GELQF'
          lqoptnb      = ilaenv( 1, subnam, ' ', nbi, ni, -1, -1 )
          factoptnb    = max(qroptnb, lqoptnb)
          IF( algo.EQ.'TRD' ) THEN
             IF( stag.EQ.'2STAG' ) THEN
                lwork = ni*nbi + ni*max(nbi+1,factoptnb)+ max(2*nbi*nbi, nbi*nthreads)+ (nbi+1)*ni
             ELSE IF( (stag.EQ.'HE2HB').OR.(stag.EQ.'SY2SB') ) THEN
                lwork = ni*nbi + ni*max(nbi,factoptnb) + 2*nbi*nbi
             ELSE IF( (stag.EQ.'HB2ST').OR.(stag.EQ.'SB2ST') ) THEN
                lwork = (2*nbi+1)*ni + nbi*nthreads
             ENDIF
          ELSE IF( algo.EQ.'BRD' ) THEN
             IF( stag.EQ.'2STAG' ) THEN
                lwork = 2*ni*nbi + ni*max(nbi+1,factoptnb)+ max(2*nbi*nbi, nbi*nthreads)+ (nbi+1)*ni
             ELSE IF( stag.EQ.'GE2GB' ) THEN
                lwork = ni*nbi + ni*max(nbi,factoptnb) + 2*nbi*nbi
             ELSE IF( stag.EQ.'GB2BD' ) THEN
                lwork = (3*nbi+1)*ni + nbi*nthreads
             ENDIF
          ENDIF
          lwork = max( 1, lwork )

          IF( lwork.GT.0 ) THEN
             iparam2stage = lwork
          ELSE
             iparam2stage = -1
          ENDIF

       ELSE IF ( ispec .EQ. 21 ) THEN

          iparam2stage = nxi
       ENDIF

       END FUNCTION iparam2stage

       SUBROUTINE dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

       DOUBLE PRECISION ALPHA,BETA
       INTEGER K,LDA,LDB,LDC,M,N
       CHARACTER TRANSA,TRANSB

       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)

       !LOGICAL LSAME

       !EXTERNAL xerbla

       INTRINSIC max

       DOUBLE PRECISION TEMP
       INTEGER I,INFO,J,L,NROWA,NROWB
       LOGICAL NOTA,NOTB

       DOUBLE PRECISION ONE,ZERO
       parameter(one=1.0d+0,zero=0.0d+0)

       nota = lsame(transa,'N')
       notb = lsame(transb,'N')
       IF (nota) THEN
           nrowa = m
       ELSE
           nrowa = k
       END IF
       IF (notb) THEN
           nrowb = k
       ELSE
           nrowb = n
       END IF

       info = 0
       IF ((.NOT.nota) .AND. (.NOT.lsame(transa,'C')) .AND. + (.NOT.lsame(transa,'T'))) THEN
           info = 1
       ELSE IF ((.NOT.notb) .AND. (.NOT.lsame(transb,'C')) .AND. + (.NOT.lsame(transb,'T'))) THEN
           info = 2
       ELSE IF (m.LT.0) THEN
           info = 3
       ELSE IF (n.LT.0) THEN
           info = 4
       ELSE IF (k.LT.0) THEN
           info = 5
       ELSE IF (lda.LT.max(1,nrowa)) THEN
           info = 8
       ELSE IF (ldb.LT.max(1,nrowb)) THEN
           info = 10
       ELSE IF (ldc.LT.max(1,m)) THEN
           info = 13
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('DGEMM ',info)
           RETURN
       END IF

       IF ((m.EQ.0) .OR. (n.EQ.0) .OR. +    (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN

       IF (alpha.EQ.zero) THEN
           IF (beta.EQ.zero) THEN
               DO 20 j = 1,n
                   DO 10 i = 1,m
                       c(i,j) = zero
    10             CONTINUE
    20         CONTINUE
           ELSE
               DO 40 j = 1,n
                   DO 30 i = 1,m
                       c(i,j) = beta*c(i,j)
    30             CONTINUE
    40         CONTINUE
           END IF
           RETURN
       END IF

       IF (notb) THEN
           IF (nota) THEN

               DO 90 j = 1,n
                   IF (beta.EQ.zero) THEN
                       DO 50 i = 1,m
                           c(i,j) = zero
    50                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 60 i = 1,m
                           c(i,j) = beta*c(i,j)
    60                 CONTINUE
                   END IF
                   DO 80 l = 1,k
                       temp = alpha*b(l,j)
                       DO 70 i = 1,m
                           c(i,j) = c(i,j) + temp*a(i,l)
    70                 CONTINUE
    80             CONTINUE
    90         CONTINUE
           ELSE

               DO 120 j = 1,n
                   DO 110 i = 1,m
                       temp = zero
                       DO 100 l = 1,k
                           temp = temp + a(l,i)*b(l,j)
   100                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   110             CONTINUE
   120         CONTINUE
           END IF
       ELSE
           IF (nota) THEN

               DO 170 j = 1,n
                   IF (beta.EQ.zero) THEN
                       DO 130 i = 1,m
                           c(i,j) = zero
   130                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 140 i = 1,m
                           c(i,j) = beta*c(i,j)
   140                 CONTINUE
                   END IF
                   DO 160 l = 1,k
                       temp = alpha*b(j,l)
                       DO 150 i = 1,m
                           c(i,j) = c(i,j) + temp*a(i,l)
   150                 CONTINUE
   160             CONTINUE
   170         CONTINUE
           ELSE

               DO 200 j = 1,n
                   DO 190 i = 1,m
                       temp = zero
                       DO 180 l = 1,k
                           temp = temp + a(l,i)*b(j,l)
   180                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   190             CONTINUE
   200         CONTINUE
           END IF
       END IF

       RETURN

       END SUBROUTINE dgemm

       SUBROUTINE xerbla( SRNAME, INFO )

       CHARACTER*(*)      SRNAME
       INTEGER            INFO

       INTRINSIC          len_trim

       WRITE( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info

       stop

  9999 FORMAT( ' ** On entry to ', a, ' parameter number ', i2, ' had ', $  'an illegal value' )

       END SUBROUTINE xerbla

       LOGICAL FUNCTION lsame(CA,CB)

       CHARACTER ca,cb

       INTRINSIC ichar

       INTEGER inta,intb,zcode

       lsame = ca .EQ. cb
       IF (lsame) RETURN

       zcode = ichar('Z')

       inta = ichar(ca)
       intb = ichar(cb)

       IF (zcode.EQ.90 .OR. zcode.EQ.122) THEN

           IF (inta.GE.97 .AND. inta.LE.122) inta = inta - 32
           IF (intb.GE.97 .AND. intb.LE.122) intb = intb - 32

       ELSE IF (zcode.EQ.233 .OR. zcode.EQ.169) THEN

           IF (inta.GE.129 .AND. inta.LE.137 .OR. + inta.GE.145 .AND. inta.LE.153 .OR. + inta.GE.162 .AND. inta.LE.169) inta = inta + 64
           IF (intb.GE.129 .AND. intb.LE.137 .OR. + intb.GE.145 .AND. intb.LE.153 .OR. + intb.GE.162 .AND. intb.LE.169) intb = intb + 64

       ELSE IF (zcode.EQ.218 .OR. zcode.EQ.250) THEN

           IF (inta.GE.225 .AND. inta.LE.250) inta = inta - 32
           IF (intb.GE.225 .AND. intb.LE.250) intb = intb - 32
       END IF
       lsame = inta .EQ. intb

       END FUNCTION lsame

       end module LinearAlgebraLowLevel