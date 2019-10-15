c     vsn 1.0  B91012: cloned from fbt.f
c     vsn 1.0  B91014: added mdetID-unwise_objid xref output
c
      character*500 InFNam, OutFnam, XrefNam, NumStr
      character*16, allocatable :: objid(:)
      character*99  TmpStr
      real*4        flux(2), dflux(2), SNRfac, zeroes(4), StDev, avg,
     +              dfarray(21) 
      real*8        ra, dec, nullval, sum, sumsq
      integer*4     nargs, iargc, nOut, FileID, status,
     +              readwrite, blocksize, nRows, nCols, nRA, nDec,
     +              nflux, ndflux, hdutype, k, felem, nelems, stat(5),
     +              narg, nobjid
      integer       n, n0, n1, n2
      logical*4     dbg, doXref
      logical       anynull
      integer*4,    allocatable :: ndex(:)
      real*4,       allocatable :: fluxes(:), dfluxes(:)
      real*8,       allocatable :: RAs(:), Decs(:)
c      
      data nOut/0/, dbg/.false./, narg/2/, zeroes/4*0.0/,
     +     doXref/.false./ 
c-----------------------------------------------------------------------
c
      nargs = iargc()
      if (nargs .lt. 2) then
        print *,'fbt2mdet vsn 1.0  B91014'
        print *,'usage: fbt2mdet infile outfile <xref|dbg>'
        print *
        print *,
     +    'where: infile  is an unWISE FITS binary table file'
        print *,
     +    '               containing a crowdsource source list'
        print *,'       outfile is the corresponding mdet list'
        print *,'       xref    is a cross-reference table for'
        print *,'               mdetID and unwise_objid'
        print *
        stop
      end if
      call getarg(1,InFNam)
      if (Access(InFNam(1:lnblnk(InFNam)),' ') .ne. 0) then
        print *,'File not found: ',InFNam(1:lnblnk(InFNam))
        call exit(64)
      end if
      call getarg(2,OutFnam)
10    if (nargs .gt. narg) then
        narg = narg + 1
        call getarg(narg, NumStr)
        if ((NumStr .eq. 'dbg') .or. (NumStr .eq. '-d')) then
          dbg = .true.
          print *,'Debug mode enabled'
        else
          XrefNam = NumStr
          doXref = .true.
        end if
        go to 10
      end if
c
c-----------------------------------------------------------------------
c
      status = 0
      call ftgiou(FileID,status)         ! get unit number
      readwrite = 0                      ! open FITS file
      call ftopen(FileID,InFNam,readwrite,blocksize,status)
      if (status .ne. 0) then
        write(6,'(a)') 'ERROR: Could not open '//trim(InFNam)
        call exit(64)
      endif
      call ftmahd(FileID,2,hdutype,status) ! move to header #2
      if (dbg) print *,'hdutype for nhdu-2, status:     ',
     +                  hdutype, status
      status = 0
      call ftgnrw(FileID,nRows,status)     ! get #rows
      if (status .ne. 0) go to 3001
      if (nRows .lt. 21) then
        print *,'ERROR: not enough data rows'
        print *,'       nRows =', nRows
        call exit(64)
      end if
      if (dbg) print *,'No. of rows in table and status:', nRows, status
c
      call ftgncl(FileID,nCols,status)     ! get #cols
      if (status .ne. 0) go to 3002
      if (dbg) print *,'No. of cols in table and status:', nCols, status
c
      call ftgcno(FileID, .false., 'ra      ', nRA,  status)  ! col# for RA
      if (status .ne. 0) go to 3003
      if (dbg) print *,'Col. no and status for ra:      ', nRA, status
c
      call ftgcno(FileID, .false., 'dec     ', nDec, status)  ! col# for Dec
      if (status .ne. 0) go to 3004
      if (dbg) print *,'Col. no and status for dec:     ', nDec, status
c
      call ftgcno(FileID, .false., 'flux    ', nflux, status) ! col# for flux
      if (status .ne. 0) go to 3005
      if (dbg) print *,'Col. no and status for flux:    ', nflux,status
c
      call ftgcno(FileID, .false., 'dflux   ', ndflux, status)! col# for dflux
      if (status .ne. 0) go to 3006
      if (dbg) print *,'Col. no and status for dflux:   ', ndflux,status
c
      if (doXref) then
        call ftgcno(FileID, .false., 'unwise_objid', nobjid, status)! col# for unwise_objid
        if (status .ne. 0) go to 3006
        if (dbg) print *,'Col. no and status for unwise_objid:',
     +           nobjid, status
      end if
c
      allocate(RAs(nRows))
      if (.not.allocated(RAs)) then
        print *,'ERROR: allocation of RAs failed'
        print *,'       no. elements =',nRows
        call exit(64)
      end if
c
      allocate(Decs(nRows))
      if (.not.allocated(Decs)) then
        print *,'ERROR: allocation of Decs failed'
        print *,'       no. elements =',nRows
        call exit(64)
      end if
c
      allocate(fluxes(nRows))
      if (.not.allocated(fluxes)) then
        print *,'ERROR: allocation of fluxes failed'
        print *,'       no. elements =',nRows
        call exit(64)
      end if
c
      allocate(dfluxes(nRows))
      if (.not.allocated(dfluxes)) then
        print *,'ERROR: allocation of dfluxes failed'
        print *,'       no. elements =',nRows
        call exit(64)
      end if
c
      allocate(ndex(nRows))
      if (.not.allocated(ndex)) then
        print *,'ERROR: allocation of ndex failed'
        print *,'       no. elements =',nRows
        call exit(64)
      end if
c
      if (doXref) then
        allocate(objid(nRows))
        if (.not.allocated(objid)) then
          print *,'ERROR: allocation of objid failed'
          print *,'       no. elements =',nRows
          call exit(64)
        end if
      end if
c          
      do 100 n = 1, nRows
c      
        felem  = 1
        nelems = 1
        k = 1
        call ftgcvd(FileID,nRA,n,felem,nelems,0.0d0,RA,anynull,status)
        if (status .ne. 0) go to 3000
        stat(1) = status
        RAs(n) = RA
        k = 2
        call ftgcvd(FileID,nDec,n,felem,nelems,0.0d0,Dec,anynull,status)
        if (status .ne. 0) go to 3000
        stat(2) = status
        Decs(n) = Dec
        if (doXref) then
          k = 5
          call FTGCvS(FileID,nobjid,n,felem,nelems,' ',TmpStr,
     +                anynull,status)
          if (status .ne. 0) go to 3000
          stat(5) = status
          objid(n) = AdjustL(TmpStr)
        end if
        nelems = 2
        k = 3
        call ftgcve(FileID,nflux,n,felem,nelems,0.0d0,flux,
     +              anynull,status)
        if (status .ne. 0) go to 3000
        stat(3) = status
        fluxes(n) = sqrt(flux(1)**2 + flux(2)**2)
        k = 4
        call ftgcve(FileID,ndflux,n,felem,nelems,0.0d0,dflux,
     +              anynull,status)
        if (status .ne. 0) go to 3000
        stat(4) = status
        dfluxes(n) = sqrt(dflux(1)**2 + dflux(2)**2)
c        
        if (dbg .and. ((n .le. 10) .or. (n .ge. nRows-10))) then
          print *,'n,felem,ra,status:   ', n, felem, RA,  stat(1)
          print *,'n,felem,dec,status:  ', n, felem, Dec, stat(2)
          print *,'n,felem,flux,status: ', n, felem, flux,  stat(3)
          print *,'n,felem,dflux,status:', n, felem, dflux, stat(4)
          if (doXref) 
     +    print *,'n,felem,objid,status:', n, felem, objid(n), stat(5)
          print *,'n,fluxes: ',fluxes(n)
          print *,'n,dfluxes:',dfluxes(n)
          print *
        end if
100   continue      
c     
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
c
      call TJISORT(nRows,fluxes,ndex)
      sum   = 0.0d0
      sumsq = 0.0d0
      n1 = nRows/2 - 10
      n2 = nRows/2 + 10
      n0 = 0
      do 120 n = n1, n2
        n0 = n0 + 1
        dfarray(n0) = dfluxes(ndex(n))
        sum   = sum   + dfluxes(ndex(n))
        sumsq = sumsq + dfluxes(ndex(n))**2
        if (dbg) print *,'n, ndex(n), fluxes(ndex(n)) dfluxes(ndex(n)):',
     +                    n, ndex(n), fluxes(ndex(n)), dfluxes(ndex(n))
120   continue
c
      call TJSORT(21,dfarray)
      avg = sum/21.0d0
      SNRfac = 1.0/dfarray(11)
      StDev  = dsqrt(dabs(sumsq/21.0d0 - avg**2))
      print *,'Median Flux:  ', fluxes(ndex(nRows/2))
      print *,'Average dFlux:', avg
      print *,'StdDev(dFlux):', StDev
      print *,'Median dFlux: ', dfarray(11)
      print *,'SNRfac:       ', SNRfac
      print *,'Median SNR:   ', SNRfac*fluxes(ndex(nRows/2))
c
      open (12, file = OutFnam)
      write(12,'(''\Nsrcs =   '',i7)') nRows
      write(12,'(a)') '| Src  |    RA    |    Dec   |   SNR   |'
     +              //' Rsat1 | Rsat2 | Rsat3 | Rsat4 |'
      write(12,'(a)') '|  i   |    r     |     r    |    r    |'
     +              //'   r   |   r   |   r   |   r   |'
      write(12,'(a)') '|      |   deg    |    deg   |         |'
     +              //' arcsec| arcsec| arcsec| arcsec|'
c
      if (doXref) then
        open (14, file = XrefNam)
        write(14,'(''\Nsrcs =   '',i7)') nRows
        write(14,'(a)') '|mdetID|  unwise_objid  |'
        write(14,'(a)') '|  int |      char      |'
      end if
c
      do 1000 n0 = 1, nRows
        n = nRows-n0+1
        write(12,'(i7,2f11.5,1pe10.3,4f8.2)')
     +             n0, RAs(ndex(n)), Decs(ndex(n)),
     +             SNRfac*fluxes(ndex(n)), zeroes
        nOut = nOut + 1
        if (doXref) write (14,'(i7,1x,a16)') n0, objid(ndex(n))
1000  continue
c
      print *,'No. of rows read:   ', nRows
      print *,'No. of rows written:', nOut
c
      stop
c
3000  print *,'ERROR reading table row no.', n,
     +        '; parameter', k, '; status =',status
      call exit(64)
3001  print *,'ERROR: unable to get @rows; status =',status
      call exit(64)
3002  print *,'ERROR: unable to get #cols; status =',status
      call exit(64)
3003  print *,'ERROR: ra not found; status =',status
      call exit(64)
3004  print *,'ERROR: dec not found; status =',status
      call exit(64)
3005  print *,'ERROR: flux not found; status =',status
      call exit(64)
3006  print *,'ERROR: dflux not found; status =',status
      call exit(64)
      stop
c      
      end
c      
c=======================================================================
c                                  Index sort for real*4 array
c                                  from Numerical Recipes via T. Jarrett
      SUBROUTINE TJISORT(N,RA,ND)
c
      Integer*4 N,L,IR,J,I,ND(N),NDA
      Real*4 RA(N),RRA
c
      if (n .lt. 1) return
      Do 5 I = 1, N
        ND(I) = I
5     Continue
      if (n .lt. 2) return
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(ND(L))
          NDA = ND(L)
        ELSE
          RRA=RA(ND(IR))
          NDA = ND(IR)
          ND(IR)=ND(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            ND(1)=NDA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(ND(J)).LT.RA(ND(J+1)))J=J+1
          ENDIF
          IF(RRA.LT.RA(ND(J)))THEN
            ND(I)=ND(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        ND(I) = NDA
      GO TO 10
      END
c
C=======================================================================
c                                  Sort for real*4 array
c                                  from Numerical Recipes via T. Jarrett
      SUBROUTINE TJSORT(N,RA)
c
      Integer*4 N,L,IR,J,I
      Real*4 RA(N),RRA
c
      if (n .lt. 2) return
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
