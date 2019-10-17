c     vsn 0.9  B91012: cloned from fbt.f
c     vsn 1.0  B91014: added mdetID-unwise_objid xref output
c     vsn 1.0  B91016: changed SNR computation per-band avg dflux)
c
      character*500 InFNam, OutFnam, XrefNam, NumStr
      character*16, allocatable :: objid(:)
      character*99  TmpStr
      real*4        flux(2), dflux(2), SNRfac1, SNRfac2, zeroes(4),
     +              StDev, avg, AvgFrac
      real*8        ra, dec, nullval, sum, sumsq
      integer*4     nargs, iargc, nOut, FileID, status,
     +              readwrite, blocksize, nRows, nCols, nRA, nDec,
     +              nflux, ndflux, hdutype, k, felem, nelems, stat(5),
     +              narg, nobjid, k1, k2
      integer       n, n0, n1, n2, nMid
      logical*4     dbg, doXref
      logical       anynull
      integer*4,    allocatable :: ndex(:)
      real*4,       allocatable :: SNRs(:), flux1s(:), dflux1s(:),
     +                                      flux2s(:), dflux2s(:)
      real*8,       allocatable :: RAs(:), Decs(:)
c      
      data nOut/0/, dbg/.false./, narg/2/, zeroes/4*0.0/,
     +     doXref/.false./, AvgFrac/0.005/ 
c-----------------------------------------------------------------------
c
      nargs = iargc()
      if (nargs .lt. 2) then
        print *,'fbt2mdet vsn 1.0  B91016'
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
      allocate(SNRs(nRows))
      if (.not.allocated(SNRs)) then
        print *,'ERROR: allocation of SNRs failed'
        print *,'       no. elements =',nRows
        call exit(64)
      end if
c
      allocate(flux1s(nRows))
      if (.not.allocated(flux1s)) then
        print *,'ERROR: allocation of flux1s failed'
        print *,'       no. elements =',nRows
        call exit(64)
      end if
c
      allocate(dflux1s(nRows))
      if (.not.allocated(dflux1s)) then
        print *,'ERROR: allocation of dflux1s failed'
        print *,'       no. elements =',nRows
        call exit(64)
      end if
c
      allocate(flux2s(nRows))
      if (.not.allocated(flux2s)) then
        print *,'ERROR: allocation of flux2s failed'
        print *,'       no. elements =',nRows
        call exit(64)
      end if
c
      allocate(dflux2s(nRows))
      if (.not.allocated(dflux2s)) then
        print *,'ERROR: allocation of dflux2s failed'
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
      n1 = 0
      n2 = 0
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
        flux1s(n)  = flux(1)
        flux2s(n)  = flux(2)
        if (flux(1) .gt. 0.0) n1 = n1 +1
        if (flux(2) .gt. 0.0) n2 = n2 +1
        k = 4
        call ftgcve(FileID,ndflux,n,felem,nelems,0.0d0,dflux,
     +              anynull,status)
        if (status .ne. 0) go to 3000
        stat(4) = status
        dflux1s(n) = dflux(1)
        dflux2s(n) = dflux(2)
c        
        if (dbg .and. ((n .le. 5) .or. (n .ge. nRows-5))) then
          print *,'n,felem,ra,status:   ', n, felem, RA,  stat(1)
          print *,'n,felem,dec,status:  ', n, felem, Dec, stat(2)
          print *,'n,felem,flux,status: ', n, felem, flux,  stat(3)
          print *,'n,felem,dflux,status:', n, felem, dflux, stat(4)
          if (doXref) 
     +    print *,'n,felem,objid,status:', n, felem, objid(n), stat(5)
          print *,'n,fluxes: ' ,flux1s(n), flux2s(n)
          print *,'n,dfluxes:',dflux1s(n), dflux2s(n)
          print *
        end if
100   continue      
c     
      call ftclos(FileID, status)
      call ftfiou(FileID, status)
c
      print *,'No. of W1 detections:', n1
      print *,'No. of W2 detections:', n2
c
      call TJISORT(nRows,flux1s,ndex)
      nMid = (2*nRows-n1)/2            ! midpoint of nonzero W1 range
      k1 = nMid - NInt(AvgFrac*float(n1))
      k2 = nMid + NInt(AvgFrac*float(n1))
      if (k1 .le. nRows-n1) k1 = nRows-n1 + 1
      if (k2 .gt. nRows)    k2 = nRows
      sum   = 0.0d0
      sumsq = 0.0d0
      do 120 n = k1, k2
        sum   = sum   + dflux1s(ndex(n))
        sumsq = sumsq + dflux1s(ndex(n))**2
120   continue
      avg = sum/float(k2-k1+1)
      SNRfac1 = 1.0/avg
      StDev  = dsqrt(dabs(sumsq/float(k2-k1+1) - avg**2))
      print *
      print *,'Median W1 Flux:  ', flux1s(ndex(nMid))
      print *,'Average W1 dFlux:', avg
      print *,'StdDev(W1 dFlux):', StDev
      print *,'Median W1 dFlux: ', dflux1s(ndex(nMid))
      print *,'SNRfac1:         ', SNRfac1
      print *,'Median W1 SNR:   ', SNRfac1*flux1s(ndex(nMid))
c
      call TJISORT(nRows,flux2s,ndex)
      nMid = (2*nRows-n2)/2            ! midpoint of nonzero W1 range
      k1 = nMid - NInt(AvgFrac*float(n2))
      k2 = nMid + NInt(AvgFrac*float(n2))
      if (k1 .le. nRows-n2) k1 = nRows-n2 + 1
      if (k2 .gt. nRows)    k2 = nRows
      sum   = 0.0d0
      sumsq = 0.0d0
      do 140 n = k1, k2
        sum   = sum   + dflux2s(ndex(n))
        sumsq = sumsq + dflux2s(ndex(n))**2
140   continue
      avg = sum/float(k2-k1+1)
      SNRfac2 = 1.0/avg
      StDev  = dsqrt(dabs(sumsq/float(k2-k1+1) - avg**2))
      print *
      print *,'Median W2 Flux:  ', flux2s(ndex(nMid))
      print *,'Average W2 dFlux:', avg
      print *,'StdDev(W2 dFlux):', StDev
      print *,'Median W2 dFlux: ', dflux2s(ndex(nMid))
      print *,'SNRfac2:         ', SNRfac2
      print *,'Median W2 SNR:   ', SNRfac2*flux1s(ndex(nMid))
c
      SNRs = 0.0
      do 200 n = 1, nRows
        SNRs(n) = sqrt((SNRFac1*flux1s(n))**2
     +                +(SNRFac2*flux2s(n))**2)
200   continue
      call TJIsort(nRows,SNRs,ndex)
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
     +             SNRs(ndex(n)), zeroes
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
