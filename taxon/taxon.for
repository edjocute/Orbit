************************************************************************
c      PROGRAM taxon
      SUBROUTINE taxon(t,x,v,np,jsub,jdim,jcla,jcl,jpan,jlin,jcom,
     \                  arch,code,peri)
*
* Orbit classifier.
*
* Version: 16.2 (24/03/2014)
*
* This program is written in FORTRAN-77 with the non-standard 
* extensions ENDDO, DOWHILE and INCLUDE.
* The program requires the companion file "taxon.def", which contains 
* parameters. The default values of those parameters are recommended.
* Main hypothesis: the major axis of the potential (or any equivalent
* measure of it) is on the x-axis (first coordinate).
*
* USAGE
*
* I) INPUT
*
*   'taxon' may be used as a PROGRAM or as a SUBROUTINE, the only
*   difference being how to obtain the few input parameters:
*   a) As a PROGRAM: the parameters are read from the file 'taxon.par',
*      which should have two records:
*      First: an alphanumeric constant, with the name of the file
*         containing the orbit. This name, with its extension (if any)
*         removed, will also be the prefix of all output files (with
*         different extensions in each case).
*      Second: five integers: 
*         1st: =2 indicates a 2D potential; = 3 a 3D potential
*         2nd: =1 dump results in file (.cla); = 0 do not.
*         3rd: =1 simplified dump of results in file (.cl); = 0 do not.
*         4th: =1 dump messages to the screen; = 0 do not.
*         5th: =1 dump line spectra to file (.lix, .llx,...); 0= do not.
*
*   b) As a SUBROUTINE: the same parameters as before are provided by 
*      the calling progran through the COMMON /input/. The order of the
*      variables is as follows:
*
*      COMMON /input/jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch
*
*      where
*         INTEGER jsub       ! Signals taxon is being called as a 
*                              subroutine; please set it non-zero.
*         INTEGER jdim       ! 2D or 3D potential.
*         INTEGER jcla       ! Dump results to file
*         INTEGER jcl        ! Dump simplified results to file
*         INTEGER jpan       ! Dump messages to the screen
*         INTEGER jlin       ! Dump line spectra to file
*         INTEGER jcom       ! The orbit is given through a COMMON
*         CHARACTER arch*30  ! The orbit is read from the file arch
*
*      If jcom <> 0, the orbit is given through a COMMON instead of 
*      reading it from a file. This COMMON must be:
*      COMMON /orbita/t(0:nn-1),x(3,0:nn-1),v(3,0:nn-1),np 
*      provided by the calling program, where t, x and v are DOUBLE 
*      PRECISION arrays, and np is an INTEGER containing the actual 
*      number of points. (Note that the dimension of these arrays is nn,
*      as declared in taxon.def, not np.) 
*
*   c) The file containing the orbit: this file is user-defined. The 
*      only requirements are that it must contain times, positions and
*      velocities of the particle at equidistant times, and that the 
*      number of positions must be a power of 2.
*      Note that the z coordinates and velocities are not needed if 
*      jdim=2. The last routine in this program ("orbitin") should be 
*      modified according to the format of the input file. If taxon is
*      used as a subroutine, this file may be ignored completely (see 
*      item b)). 
*
* II) OUTPUT
*
*   The output is driven by the input parameters described above.
*
*   a) The results of the classification (switch 'jcla') are driven by
*      the subroutine "salida", which may be modified as desired to
*      fulfill requirements of the user. The output (file with extension
*      .cla) is plenty of banners, so it is auto-explained.
*
*   b) An alternative simplified output (switch 'jcl') is output to a
*      file with extension .cl. By default, this option do not destroy
*      previous records of the .cl file. Each record of this file will
*      contain:
*         classification code (see below), integers of the resonances
*         (if any).
*
*   c) Messages of progress in the computation, warnings, and results,
*      are dumped to the screen (switch 'jpan').
*
*   d) When used as a subroutine, the results may also be obtained by
*      putting the following COMMON /results/ in the calling program:
*
*      COMMON /results/ffund(4),code,lic1(3),lic2(3),licf(3),peri
*
*      where:
*      DOUBLE PRECISION ffund(4) are the main frequencies;
*        0 if any is absent
*      INTEGER code is the classification code (see below). 
*      INTEGER lic1(3) are the integers of a first resonance.
*      INTEGER lic2(3) are the integers of a second resonance.
*      INTEGER licf(3) are the integers of the full resonance (when 2
*              resonances are present)
*      INTEGER peri contains the number of orbital periods.
*
*   e) Line spectra (switch 'jlin'): a file is generated for each 
*      coordinate (with extensions .lix, .liy and .liz), with the 
*      following format:
*
*      record 1: frequency  0          0      [of line 1]
*      record 2: frequency  amplitude  phase  [of line 1]
*      record 3: [blank]
*      record 4: frequency  0          0      [of line 2]
*      etc.
*
*      This format allows to plot the lines with the convention that
*      a blank line should discontinue the lines of the plot. The phase
*      is given in sexagesimal degrees. Also, logarithmic line spectra
*      are dumped to files (with extensions .llx, .lly and .llz), with
*      the following format:
*
*      record 1: frequency -10              freq [rad/u.t.] [of line 1]
*      record 2: frequency LOG10(amplitude) freq [rad/u.t.] [of line 1]
*      record 3: [blank]
*      record 4: frequency -10              freq [rad/u.t.] [of line 2]
*      etc.
*
*
*   TABLE OF CODES
*
*   1st column: The codes of the final classification, which appear in
*      the .cla and .cl files, and in the COMMON /results/. They are 
*      numbers of three digits. The first one indicates the number of
*      independent (fundamental) frequencies of the orbit; the second, 
*      the number of resonances; the third, the morphology (0 = box, 
*      1 = x tube; 2 = y tube; 3 = z tube).
*      Note I: don't mix up main frequencies with fundamental
*       frequencies. The former are those with greater amplitude, with
*       which the resonances are made up. The latter are the linearly
*       independent frequencies.
*      Note II: y tubes may appear if the orientation of the potential
*       is such that the y coordinate is the minor axis.
*   
*   2nd column: "dimension" stands for the dimension of the subspace
*      into which the orbit moves (given by the number of fundamental
*      frequencies). The names in each case are as follows:
*      Dimension   Name when a 2D potential  Name when a 3D potential
*      1           closed                    closed
*      2           open                      thin
*      3           irregular                 open
*      4           -                         irregular
*
*   3rd column: the number of resonances allows to determine the family
*      to which the orbit belongs. The parent of the orbit is listed.
*
*   4th, 5th, 6th columns: morphology, number of fundamental 
*      frequencies, and number of resonances.
*
* code  dimension parent  morphol   nff  nres  Obs
*-----------------------------------------------------------------------
* 000   not classified     -        -    -
* 100   closed      -     box       1    0    1D (linear) orbit
* 110   closed      -     box       1    1    2D orbit
* 111   closed      -     loop      1    1    2D orbit
* 120   closed      -     box       1    2
* 121   closed      -     x tube    1    2
* 122   closed      -     y tube    1    2
* 123   closed      -     z tube    1    2
* 200   open        -     box       2    0    2D potential
* 210   open      closed  box       2    1    2D potential
*  or   thin        -     box       2    1
* 211   open      closed  loop      2    1    2D potential
*  or   thin        -     x tube    2    1
* 212   thin        -     y tube    2    1
* 213   thin        -     z tube    2    1
* 220   thin      closed  box       2    2
* 221   thin      closed  x tube    2    2
* 222   thin      closed  y tube    2    2
* 223   thin      closed  z tube    2    2
* 300   irregular   -      -       >2    -    2D potential
*  or   open        -     box       3    0
* 310   open      thin    box       3    1
* 311   open      thin    x tube    3    1
* 312   open      thin    y tube    3    1
* 313   open      thin    z tube    3    1
* 320   open      closed  box       3    2
* 321   open      closed  x tube    3    2
* 322   open      closed  y tube    3    2
* 323   open      closed  z tube    3    2
* 400   irregular    -     -       >3    -
*
* Please send comments, bugs, etc. to: ddc@fcaglp.unlp.edu.ar
************************************************************************
      INCLUDE 'taxon.def'
      INTEGER ju
* Orbit
      INTEGER np
      DOUBLE PRECISION t(0:nn-1),x(3,0:nn-1),v(3,0:nn-1)
*      COMMON /orbita/t,x,v,np
*      SAVE /orbita/
Cf2py intent(in) np
Cf2py intent(in) t,x,v
* Input
      INTEGER jsub,jdim,jcla,jcl,jpan,jlin,jcom
      CHARACTER arch*30
c      COMMON /input/jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch
Cf2py intent(in) jsub,jdim,jcla,jcl,jpan,jlin,jcom
Cf2py intent(in) arch
* Output
      INTEGER code,peri
Cf2py intent(out) code,peri

* Input data
      CALL entrada(ju,np,t,x,v,jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch)
* Line spectra
      CALL nesvor(np,t,x,v,jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch)
* Clean up the spectra
      CALL limpia(jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch)
* Classify the orbit
      CALL clasificador(jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch,
     \                  code,peri)
* Output
      CALL salida(ju,np,t,x,v,jsub,jdim,jcla,jcl,jpan,jlin,jcom,
     \            arch,code,peri)
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                        CLASSIFICATION SECTION                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

************************************************************************
      SUBROUTINE clasificador(jsub,jdim,jcla,jcl,jpan,jlin,jcom,
     \                        arch,code,peri)
* It classifies an orbit from its line spectra.
* Output: ffund, code, lic1, lic2, licf (COMMON /results/).
************************************************************************
      INCLUDE 'taxon.def'
* Line spectra
      INTEGER nl(3)
      DOUBLE PRECISION amp(3,maxlin),fas(3,maxlin),fre(3,maxlin)
      COMMON /especlineas/amp,fre,fas,nl
      SAVE /especlineas/
* Input
      INTEGER jsub,jdim,jcla,jcl,jpan,jlin,jcom
      CHARACTER arch*30
*      COMMON /input/jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch
* Output
      INTEGER code,lic1(3),lic2(3),licf(3),peri
      DOUBLE PRECISION ffund(4)
      COMMON /results/ffund,lic1,lic2,licf
      SAVE /results/
Cf2py intent(out) code,peri

      INTEGER l,k,ndim,nff,nres,u,nirr,i1,i2,j,ind(3*maxlin),r1,r2,r3,
     & r4,r5,r6,cfm
      DOUBLE PRECISION fnew,aord(3*maxlin),
     & ford(3*maxlin),faux(3*maxlin),fmain(3)

      IF(jpan.NE.0)WRITE(*,*)'Clasifying the orbit...'
* Safety valve (orbits which escape linearly for example)
      IF(nl(1)+nl(2)+nl(3).EQ.0)THEN
         code=0
         RETURN
      ENDIF
* Initialization
      DO k=1,3
         lic1(k)=0
         lic2(k)=0
         licf(k)=0
* Main frequencies allow to compute the resonances; fundamental
* frequencies determine the dimension of the manifold in which the
* orbit moves.
         fmain(k)=0d0 ! main frequencies
         ffund(k)=0d0 ! fundamental frequencies
      ENDDO
      ffund(4)=0d0

* Dimension of the orbit in configuration space
      ndim=0
      DO k=1,3
         IF(nl(k).NE.0)ndim=ndim+1
      ENDDO

* Main frequencies 
      CALL principales(fmain,cfm,nirr)

* Order the lines of all coordinates by amplitude (later used in 
* searching for fundamental frequencies)
      j=0
      DO k=1,3
         DO l=1,maxlin
            j=j+1
            aord(j)=-amp(k,l)
            faux(j)=fre(k,l)
         ENDDO
      ENDDO
      CALL quicki8(3*maxlin,aord,ind)
      DO k=1,3*maxlin
         ford(k)=faux(ind(k))
      ENDDO

* Case 1D
      IF(ndim.EQ.1)THEN
* No resonances in 1D
         nres=0
* Triplet of integers; main freq. is fundamental freq. in 1D
         IF(nl(1).NE.0)THEN
            lic1(1)=1
            ffund(1)=fmain(1)
         ELSEIF(nl(2).NE.0)THEN
            lic1(2)=1
            ffund(1)=fmain(2)
         ELSEIF(nl(3).NE.0)THEN
            lic1(3)=1
            ffund(1)=fmain(3)
         ENDIF
* Search for an additional independent (fundamental) frequency.
         CALL nuevafrec2(ford,ffund(1),fnew)
         IF(fnew.NE.0d0)ffund(2)=ABS(fnew)
* Number of independent (fundamental) frequencies
         nff=1
         IF(ffund(2).NE.0d0)nff=2
* If nirr=1, there was a problem with the order of the frequencies,
* so the orbit should be chaotic in this case.
         IF(nirr.NE.0.AND.nff.NE.2)STOP'Wrong frequencies.'

* Case 2D
      ELSEIF(ndim.EQ.2)THEN
* Indexes to the coordinates that are present
         IF(nl(1).EQ.0)THEN
            i1=2
            i2=3
         ELSEIF(nl(2).EQ.0)THEN
            i1=1
            i2=3
         ELSE
            i1=1
            i2=2
         ENDIF
* Resonances between main frequencies
         CALL resonancias2d(fmain(i1),fmain(i2),r1,r2,nres)
         lic1(i1)=r1
         lic1(i2)=r2
* Fundamental frequencies
         ffund(1)=fmain(i1)
         IF(nres.EQ.0)THEN
            ffund(2)=fmain(i2)
            CALL nuevafrec3(ford,ffund(1),ffund(2),fnew)
            IF(fnew.NE.0d0)ffund(3)=ABS(fnew)
         ELSE
* Find new fundamental frequencies
            CALL nuevafrec2(ford,ffund(1),fnew)
            IF(fnew.NE.0d0)THEN
               ffund(2)=ABS(fnew)
               CALL nuevafrec3(ford,ffund(1),ffund(2),fnew)
               IF(fnew.NE.0d0)ffund(3)=ABS(fnew)
            ENDIF
         ENDIF
* Number of fundamental frequencies
         nff=1
         IF(ffund(2).NE.0d0)nff=2
         IF(ffund(3).NE.0d0)nff=3
* If nirr=1, there was a problem with the order of the frequencies,
* so the orbit should be chaotic in this case.
         IF(nirr.NE.0.AND.nff.NE.3)STOP'Wrong frequencies.'

* Case 3D
      ELSEIF(ndim.EQ.3)THEN
* Resonances between main frequencies
         CALL resonancias(fmain(1),fmain(2),fmain(3),r1,r2,r3,r4,r5,r6,
     &    nres)
         lic1(1)=r1
         lic1(2)=r2
         lic1(3)=r3
         lic2(1)=r4
         lic2(2)=r5
         lic2(3)=r6
* Search for fundamental frequencies.
* If nres=0, the three main are fundamental; we look for a fourth freq
         IF(nres.EQ.0)THEN
            ffund(1)=fmain(1)
            ffund(2)=fmain(2)
            ffund(3)=fmain(3)
            CALL nuevafrec4(ford,ffund(1),ffund(2),ffund(3),fnew)
            IF(fnew.NE.0d0)ffund(4)=ABS(fnew)
* If nres=1:
         ELSEIF(nres.EQ.1)THEN
* If the resonance is among only two of the coordinates, there are
* only two fundamental frequencies so far
            IF(lic1(1)*lic1(2)*lic1(3).EQ.0)THEN
               IF(lic1(3).EQ.0)THEN ! z independent
                  ffund(1)=fmain(1) ! one member of the pair
                  ffund(2)=fmain(3) ! the non-paired
               ELSEIF(lic1(2).EQ.0)THEN
                  ffund(1)=fmain(1)
                  ffund(2)=fmain(2)
               ELSEIF(lic1(1).EQ.0)THEN
                  ffund(1)=fmain(2)
                  ffund(2)=fmain(1)
               ENDIF
               CALL nuevafrec3(ford,ffund(1),ffund(2),fnew)
               IF(fnew.NE.0d0)THEN
                  ffund(3)=ABS(fnew)
                  CALL nuevafrec4(ford,ffund(1),ffund(2),ffund(3),fnew)
                  IF(fnew.NE.0d0)ffund(4)=ABS(fnew)
               ENDIF
            ELSE
* If the resonance is among the three coordinates: only 1 fund. freq.
               ffund(1)=fmain(1)
               CALL nuevafrec2(ford,ffund(1),fnew)
               IF(fnew.NE.0d0)THEN
                  ffund(2)=ABS(fnew)
                  CALL nuevafrec3(ford,ffund(1),ffund(2),fnew)
                  IF(fnew.NE.0d0)THEN
                     ffund(3)=ABS(fnew)
                     CALL nuevafrec4(ford,ffund(1),ffund(2),ffund(3),
     &                fnew)
                     IF(fnew.NE.0d0)ffund(4)=ABS(fnew)
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF(nres.EQ.2)THEN
* If nres=2, there is only one fundamental; we look for more
            ffund(1)=fmain(1)
            CALL nuevafrec2(ford,ffund(1),fnew)
            IF(fnew.NE.0d0)THEN
               ffund(2)=ABS(fnew)
               CALL nuevafrec3(ford,ffund(1),ffund(2),fnew)
               IF(fnew.NE.0d0)THEN
                  ffund(3)=ABS(fnew)
                  CALL nuevafrec4(ford,ffund(1),ffund(2),ffund(3),fnew)
                  IF(fnew.NE.0d0)ffund(4)=ABS(fnew)
               ENDIF
            ENDIF
* Full resonance
            licf(1)=ABS(lic1(2)*lic2(3)-lic1(3)*lic2(2))
            licf(2)=ABS(lic1(1)*lic2(3)-lic1(3)*lic2(1))
            licf(3)=ABS(lic1(2)*lic2(1)-lic1(1)*lic2(2))
         ENDIF
* Number of fundamental frequencies
         nff=1
         IF(ffund(2).NE.0d0)nff=2
         IF(ffund(3).NE.0d0)nff=3
         IF(ffund(4).NE.0d0)nff=4
* If nirr=1, there was a problem with the order of the frequencies,
* so the orbit should be chaotic in this case.
*         IF(nirr.NE.0.AND.nff.NE.4)STOP'Wrong frequencies.'
      ENDIF

* Subfamilies
      IF(nres.GT.0)CALL impropia(ffund(1),cfm,lic1,lic2,licf)

* Classification code
      u=0
      code=nff*100
      IF(nres.GT.0)THEN
         code=code+nres*10
* The orbit is a resonant box unless...
         IF(jdim.EQ.2)THEN
* n:n is a loop unless the phases are the same or differ by pi
            IF(lic1(1).EQ.-lic1(2).AND.
     &       ABS(TAN(fas(1,1))-TAN(fas(2,1))).GT.1d-4)u=1
         ELSE
            IF(nres.EQ.1)THEN
            IF(lic1(2).EQ.-lic1(3))u=1
            IF(lic1(1).EQ.-lic1(3))u=2
            IF(lic1(1).EQ.-lic1(2))u=3
            ELSE
               IF(licf(2).EQ.licf(3))u=1
               IF(licf(1).EQ.licf(3))u=2
               IF(licf(1).EQ.licf(2))u=3
            ENDIF
         ENDIF
         code=code+u
      ENDIF    
      IF(nirr.NE.0.AND.nff.NE.4) code=-1
      END

************************************************************************
      SUBROUTINE principales(fmain,cfm,nirr)
      INTEGER nirr,cfm
      DOUBLE PRECISION fmain(3)
* It searches for main frequencies fmain.
* cfm: coordinate (1-3) to which fmain(1) belongs.
* nirr: if the lines are disordered, nirr=1 (irregular orbit)
************************************************************************
      INCLUDE 'taxon.def'
      INTEGER nl(3)
      DOUBLE PRECISION amp(3,maxlin),fas(3,maxlin),fre(3,maxlin)
      COMMON /especlineas/amp,fre,fas,nl
      SAVE /especlineas/

      LOGICAL is0,ltne
      INTEGER i,ind(maxlin),k,l
      DOUBLE PRECISION x1,x2,aord(maxlin),ford(3,maxlin)
* Less than and not equal
      ltne(x1,x2)=x1.LT.x2.AND..NOT.is0(x1,-x2)

      IF(metodo.EQ.1)THEN
* Use positive frequencies
         DO k=1,3
            DO l=1,maxlin
               ford(k,l)=ABS(fre(k,l))
            ENDDO
         ENDDO
      ELSE
* Compute amplitudes as if x + i vx 
         DO k=1,3
            DO l=1,maxlin
               aord(l)=-amp(k,l)*ABS(1d0-2d0*ACOS(-1d0)*fre(k,l))
            ENDDO
            CALL quicki8(maxlin,aord,ind)
            DO l=1,maxlin
               ford(k,l)=ABS(fre(k,ind(l)))
            ENDDO
         ENDDO
      ENDIF

      nirr=1 ! Irregular orbit
      IF(nl(1).GT.0)THEN
         fmain(1)=ford(1,1)
         cfm=1
      ELSEIF(nl(2).GT.0)THEN
         fmain(2)=ford(2,1)
         cfm=2
      ELSE
         fmain(3)=ford(3,1)
         cfm=3
      ENDIF

      IF(nl(2).NE.0.AND.nl(1).NE.0)THEN
         i=1
* MIN(i,nl(2)) appears because Fortran may evaluate the ltne function
* BEFORE the i.LE.nl(2)
         DOWHILE(i.LE.nl(2).AND.ltne(ford(2,1),ford(1,MIN(i,nl(2)))))
            i=i+1
         ENDDO
         IF(i.GT.nl(2))RETURN
         fmain(2)=ford(2,i)
      ENDIF
      IF(nl(2).NE.0.AND.nl(3).NE.0)THEN
         i=1
         DOWHILE(i.LE.nl(3).AND.ltne(ford(3,1),ford(2,MIN(i,nl(3)))))
            i=i+1
         ENDDO
         IF(i.GT.nl(3))RETURN
         fmain(3)=ford(3,i)
      ENDIF
      IF(nl(3).NE.0.AND.nl(1).NE.0)THEN
         i=1
         DOWHILE(i.LE.nl(3).AND.ltne(ford(3,1),ford(1,MIN(i,nl(3)))))
            i=i+1
         ENDDO
         IF(i.GT.nl(3))RETURN
         fmain(3)=ford(3,i)
      ENDIF
      nirr=0
      END

************************************************************************
      SUBROUTINE resonancias2d(f1,f2,i1,i2,nres)
      INTEGER i1,i2,nres
      DOUBLE PRECISION f1,f2
* Searching for resonances between f1 and f2
* (i1,i2) integer resonant vector; zeroes if there is not a resonance.
* nres: number of resonances (0 or 1) found.
************************************************************************
      INCLUDE 'taxon.def'
      LOGICAL is0
      DOUBLE PRECISION fi1,fi2

      nres=1
      DO i1=-1,-ifin2,-1
         fi1=i1*ABS(f1)
         DO i2=1,ifin2
            fi2=i2*ABS(f2)
* Linear combination: if zero, there is resonance
            IF(is0(fi1,fi2))RETURN
         ENDDO
      ENDDO
      nres=0
      i1=0
      i2=0
      END

************************************************************************
      SUBROUTINE resonancias(f1,f2,f3,ii1,ii2,ii3,jj1,jj2,jj3,nres)
      INTEGER ii1,ii2,ii3,jj1,jj2,jj3,nres
      DOUBLE PRECISION f1,f2,f3
* Searching for resonances between f1, f2 and f3.
* (ii1,ii2,ii3) and (jj1,jj2,jj3): integer resonant vectors; zeroes if
*    there is not resonance.
* nres: number of resonances (0, 1 or 2) found.
************************************************************************
      DOUBLE PRECISION eps
      COMMON /epsilon/eps
      SAVE /epsilon/

      INCLUDE 'taxon.def'
      INTEGER i1,i2,i3,isum,jsum,kk1,kk2,kk3,ksum,k2,k3
      DOUBLE PRECISION fi1,fi2,fi3

      nres=0
      ii1=0
      ii2=0
      ii3=0
      jj1=0
      jj2=0
      jj3=0
      DO i1=0,ifin3
         fi1=i1*f1
         DO i2=0,ifin3
         DO k2=-1,1,2
            fi2=k2*i2*f2
            DO i3=0,ifin3
            DO k3=-1,1,2
               fi3=k3*i3*f3
* Triple zero is not a solution
               IF(i1+i2+i3.NE.0)THEN
* Linear combination: if it is zero, i.e., if there is resonance
                  IF(ABS(fi1+fi2+fi3)/SQRT(fi1**2+fi2**2+fi3**2)
     &             .LT.eps)THEN
* If it is the first resonance found
                     IF(nres.EQ.0)THEN
* Computation of the three coprime integers
                        ii1=i1
                        ii2=k2*i2
                        ii3=k3*i3
                        isum=ABS(ii1)+ABS(ii2)+ABS(ii3)
                        nres=1
* If it is the second resonance found
                     ELSEIF(nres.EQ.1)THEN
* Computation of the three coprime integers
                        CALL terna(i1,k2*i2,k3*i3,jj1,jj2,jj3)
* Verification that it is not the first resonance
                        IF(.NOT.(ABS(ii1).EQ.ABS(jj1).AND.
     &                           ABS(ii2).EQ.ABS(jj2).AND.
     &                           ABS(ii3).EQ.ABS(jj3)))THEN
                           jsum=ABS(jj1)+ABS(jj2)+ABS(jj3)
                           nres=2
                        ELSE
                           jj1=0
                           jj2=0
                           jj3=0
                        ENDIF
* If two resonances were already found, we verify whether the new one
* has smaller integers
                     ELSE
* Computation of the three coprime integers
                        CALL terna(i1,k2*i2,k3*i3,kk1,kk2,kk3)
* Verification that it is not the first or the second resonance
                        IF(.NOT.(ABS(ii1).EQ.ABS(kk1).AND.
     &                           ABS(ii2).EQ.ABS(kk2).AND.
     &                           ABS(ii3).EQ.ABS(kk3)).AND.
     &                     .NOT.(ABS(jj1).EQ.ABS(kk1).AND.
     &                           ABS(jj2).EQ.ABS(kk2).AND.
     &                           ABS(jj3).EQ.ABS(kk3)))THEN
* We verify whether the integers are smaller
                           ksum=ABS(kk1)+ABS(kk2)+ABS(kk3)
                           IF(isum.GE.jsum.AND.isum.GE.ksum)THEN
                              ii1=kk1
                              ii2=kk2
                              ii3=kk3
                              isum=ksum
                           ELSEIF(jsum.GE.isum.AND.jsum.GE.ksum)THEN
                              jj1=kk1
                              jj2=kk2
                              jj3=kk3
                              jsum=ksum
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      END

************************************************************************
      SUBROUTINE nuevafrec2(ford,f1,fc)
      INCLUDE 'taxon.def'
      DOUBLE PRECISION ford(3*maxlin),f1,fc
* Search for a frequency fc (linearly) independent of f1
* ford: all lines ordered by amplitude
* fc = 0: such a frequency was not found.
************************************************************************
* Line spectra
      INTEGER nl(3)
      DOUBLE PRECISION amp(3,maxlin),fas(3,maxlin),fre(3,maxlin)
      COMMON /especlineas/amp,fre,fas,nl
      SAVE /especlineas/

      LOGICAL is0
      INTEGER i0,i1,k,k1
      DOUBLE PRECISION fi1,fi2

* Sweep lines
      DO k=1,nl(1)+nl(2)+nl(3)
         IF(ford(k).NE.f1)THEN
            fc=ford(k)
            DO i0=1,ifon
               fi1=i0*fc
* Linear combinations with +/-f1
               DO i1=0,ifun+i0*ifon/2 ! heuristic
                  DO k1=-1,1,2
                     fi2=k1*i1*f1
                     IF(is0(fi1,fi2))GOTO 1 ! new fc
                  ENDDO
               ENDDO
            ENDDO
            RETURN ! there was no l.c. --> fc indep.
         ENDIF
1        CONTINUE
      ENDDO
      fc=0d0
      END

************************************************************************
      SUBROUTINE nuevafrec3(ford,f1,f2,fc)
      INCLUDE 'taxon.def'
      DOUBLE PRECISION ford(3*maxlin),f1,f2,fc
* Search for a frequency fc (linearly) independent of f1 and f2
* ford: all lines ordered by amplitude
* fc = 0: such a frequency was not found.
************************************************************************
* Line spectra
      INTEGER nl(3)
      DOUBLE PRECISION amp(3,maxlin),fas(3,maxlin),fre(3,maxlin)
      COMMON /especlineas/amp,fre,fas,nl
      SAVE /especlineas/

      DOUBLE PRECISION eps
      COMMON /epsilon/eps
      SAVE /epsilon/

      INTEGER i0,i1,i2,k,k1,k2
      DOUBLE PRECISION fi0,fi1,fi2

* Sweep lines
      DO k=1,nl(1)+nl(2)+nl(3)
         IF(ford(k).NE.f1.AND.ford(k).NE.f2)THEN
            fc=ford(k)
            DO i0=1,ifon
               fi0=i0*fc
* Linear combinations with +/-f1 and +/-f2
               DO i1=0,ifun
                  DO k1=-1,1,2
                     fi1=k1*i1*f1
                     DO i2=0,ifun
                        DO k2=-1,1,2
                           IF(i0+ABS(i1)+ABS(i2).LE.ifun)THEN
                              fi2=k2*i2*f2
                              IF(ABS(fi0+fi1+fi2)/
     &                         SQRT(fi0**2+fi1**2+fi2**2).LT.eps)GOTO 1
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            RETURN ! there was no l.c. --> fc indep.
         ENDIF
1        CONTINUE
      ENDDO
      fc=0d0
      END

************************************************************************
      SUBROUTINE nuevafrec4(ford,f1,f2,f3,fc)
      INCLUDE 'taxon.def'
      DOUBLE PRECISION ford(3*maxlin),f1,f2,f3,fc
* Search for a frequency fc (linearly) independent of f1, f2 and f3
* ford: all lines ordered by amplitude
* fc = 0: such a frequency was not found.
************************************************************************
* Line spectra
      INTEGER nl(3)
      DOUBLE PRECISION amp(3,maxlin),fas(3,maxlin),fre(3,maxlin)
      COMMON /especlineas/amp,fre,fas,nl
      SAVE /especlineas/

      DOUBLE PRECISION eps
      COMMON /epsilon/eps
      SAVE /epsilon/

      INTEGER i0,i1,i2,i3,k,k1,k2,k3
      DOUBLE PRECISION fi0,fi1,fi2,fi3

* Sweep lines
      DO k=1,nl(1)+nl(2)+nl(3)
         IF(ford(k).NE.f1.AND.ford(k).NE.f2.AND.ford(k).NE.f3)THEN
            fc=ford(k)
            DO i0=1,ifon
               fi0=i0*fc
* Linear combinations with +/-f1, +/-f2 and +/-f3
               DO i1=0,ifun
                  DO k1=-1,1,2
                     fi1=k1*i1*f1
                     DO i2=0,ifun
                        DO k2=-1,1,2
                           fi2=k2*i2*f2
                           DO i3=0,ifun
                              DO k3=-1,1,2
                                 IF(i0+ABS(i1)+ABS(i2)+ABS(i3).LE.ifun)
     &                            THEN
                                    fi3=k3*i3*f3
                                    IF(ABS(fi0+fi1+fi2+fi3)/
     &                               SQRT(fi0**2+fi1**2+fi2**2+fi3**2)
     &                               .LT.eps)GOTO 1
                                 ENDIF
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            RETURN ! there was no l.c. --> fc indep.
         ENDIF
1        CONTINUE
      ENDDO
      fc=0d0
      END

************************************************************************
      SUBROUTINE impropia(bf1,cfm,t1,t2,tf)
      INTEGER cfm,t1(3),t2(3),tf(3)
      DOUBLE PRECISION bf1
* Search for higher rank resonances
************************************************************************
      INCLUDE 'taxon.def'
* Line spectra
      INTEGER nl(3)
      DOUBLE PRECISION amp(3,maxlin),fas(3,maxlin),fre(3,maxlin)
      COMMON /especlineas/amp,fre,fas,nl
      SAVE /especlineas/

      LOGICAL is0
      INTEGER k,l,j,uno,euclides,i,m
      DOUBLE PRECISION ampco,fm1,fm2

* cfm is the coordinate to which bf1 belongs
      k=cfm
      ampco=0.05d0*amp(k,1)

* Determine subfamily (heuristic approach): there is a subfamily if
* there is a line on the bf1's coordinate with an m:m-1 relation wrt 
* bf1, and no symmetric line wrt bf1
      DO uno=-1,1,2 ! m+1 and m-1
         DO m=2,5 ! multiplier of the resonance
* Search for a line 2:1, 3:2, etc.
            DO j=1,nl(k)
               fm1=(m+uno)*ABS(bf1)
               fm2=-m*ABS(fre(k,j))
               IF(is0(fm1,fm2).AND.amp(k,j).GT.ampco)THEN ! candidate
                  fm1=(m-uno)*ABS(bf1) ! go find a symmetric line
                  DO l=1,nl(k)
                     fm2=-m*ABS(fre(k,l))
* If there is a symmetric line with enough amplitude, then no subfamily
                     IF(is0(fm1,fm2).AND.amp(k,l).GT.ampco)GOTO 1
                  ENDDO
* There were no symmetric lines -> improper resonance
                  DO l=1,3
                     t1(l)=m*t1(l)
                     t2(l)=m*t2(l)
                     tf(l)=m*tf(l)
                  ENDDO
                  RETURN
               ENDIF
1              CONTINUE
            ENDDO
         ENDDO
      ENDDO
      ampco=0.2d0*amp(k,1)
* Determine subfamily (heuristic approach): there is a subfamily if
* there is a line on the bf1's coordinate with an (1...m-1:m) relation 
* wrt bf1, and no symmetric line wrt bf1
      DO m=9,2,-1 ! multiplier of the resonance (catch high res. first)
         DO i=1,m-1
            IF(euclides(i,m).EQ.1)THEN
               DO j=1,nl(k)
                  fm1=i*ABS(bf1)
                  fm2=m*ABS(fre(k,j))
                  IF(is0(fm1,-fm2).AND.amp(k,j).GT.ampco)THEN
                     fm1=2*ABS(bf1) ! go find a symmetric line
                     DO l=1,nl(k)
                        fm2=ABS(fre(k,l))+ABS(fre(k,j))
* If there is a symmetric line with enough amplitude, then no subfamily
                        IF(is0(fm1,-fm2).AND.amp(k,l).GT.ampco)GOTO 2
                     ENDDO
* There were no symmetric lines -> improper resonance
                     DO l=1,3
                        t1(l)=m*t1(l)
                        t2(l)=m*t2(l)
                        tf(l)=m*tf(l)
                     ENDDO
                     RETURN
                  ENDIF
2                 CONTINUE
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      END

***********************************************************************
      SUBROUTINE tubox(io,np,t,x,v)
      CHARACTER io*1
* Classifies an x-tube orbit as inner or outer, from a morphological 
* point of view. 
************************************************************************
      INCLUDE 'taxon.def'
* Orbit
      INTEGER np
      DOUBLE PRECISION t(0:nn-1),x(3,0:nn-1),v(3,0:nn-1)
*      COMMON /orbita/t,x,v,np
*      SAVE /orbita/

      INTEGER i
      DOUBLE PRECISION fb,fc,xmax,ymbor,ymcen,za

* Maximum x when z=0
      xmax=-1d0
      za=x(3,0)
      DO i=1,np-1
         IF(x(3,i)*za.LT.0d0)THEN
            xmax=MAX(xmax,ABS(x(1,i)))
         ENDIF
         za=x(3,i)
      ENDDO
* Central and border strips
      fc=xmax/5d0
      fb=xmax-fc
* Central and border ymax's
      ymcen=-1d0
      ymbor=-1d0
      za=x(3,0)
      DO i=1,np-1
         IF(x(3,i)*za.LT.0d0)THEN
            IF(ABS(x(1,i)).LT.fc)THEN
               ymcen=MAX(ymcen,ABS(x(2,i)))
            ENDIF
            IF(ABS(x(1,i)).GT.fb)THEN
               ymbor=MAX(ymbor,ABS(x(2,i)))
            ENDIF
         ENDIF
         za=x(3,i)
      ENDDO
      IF(ymcen.LT.ymbor)THEN
         io='i' ! inner x-tube
      ELSE
         io='o' ! outer x-tube
      ENDIF
      END

************************************************************************
      SUBROUTINE terna(i1,i2,i3,j1,j2,j3)
      INTEGER i1,i2,i3,j1,j2,j3
* Search for the coprime version j1,j2,j3 of the triplet i1,i2,i3.
************************************************************************
      INTEGER euclides,ia1,ia2,ia3,mcd,mcd2

      ia1=ABS(i1)
      ia2=ABS(i2)
      ia3=ABS(i3)
* If three zeroes, nothing is done
      IF(ia1+ia2+ia3.EQ.0)RETURN
      mcd=euclides(ia1,ia2)
      mcd2=euclides(ia3,mcd)
      j1=SIGN(ia1/mcd2,i1)
      j2=SIGN(ia2/mcd2,i2)
      j3=SIGN(ia3/mcd2,i3)
      END

************************************************************************
      FUNCTION euclides(u0,v0)
      INTEGER euclides,u0,v0
* Find the g.c.d between u0 and v0, using Euclids' algorithm
************************************************************************
      INTEGER aux,u,v

      u=u0
      v=v0
      DOWHILE(v.NE.0)
         aux=MOD(u,v)
         u=v
         v=aux
      ENDDO
      euclides=u
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                          EXTRACTION SECTION                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

************************************************************************
      SUBROUTINE limpia(jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch)
* Clean up the spectra by discarding frequency zero and very small
* amplitudes.
************************************************************************
      INCLUDE 'taxon.def'

* Input
      INTEGER jsub,jdim,jcla,jcl,jpan,jlin,jcom
      CHARACTER arch*30
*      COMMON /input/jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch

      INTEGER nl(3)
      DOUBLE PRECISION amp(3,maxlin),fas(3,maxlin),fre(3,maxlin)
      COMMON /especlineas/amp,fre,fas,nl
      SAVE /especlineas/

      LOGICAL is0
      INTEGER k,j
      DOUBLE PRECISION ampco

* Remove zero frequencies and small amplitude lines
      ampco=frac*MAX(amp(1,1),amp(2,1),amp(3,1))
      DO k=1,jdim
         nl(k)=0
         DO j=1,maxlin
            IF(.NOT.is0(fre(k,j),0d0).AND.amp(k,j).GT.ampco)THEN
               nl(k)=nl(k)+1
               fre(k,nl(k))=fre(k,j)
               amp(k,nl(k))=amp(k,j)
               fas(k,nl(k))=fas(k,j)
            ENDIF
         ENDDO
      ENDDO
      END

************************************************************************
      SUBROUTINE nesvor(np,t,x,v,jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch)
* It computes frequencies, amplitudes and phases from positions and
* velocities (Sidlichovsky and Nesvorny 1997, Cel. Mech. 65, 137). 
************************************************************************
      INCLUDE 'taxon.def'
* Input
      INTEGER jsub,jdim,jcla,jcl,jpan,jlin,jcom
      CHARACTER arch*30
*      COMMON /input/jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch
* Line spectra
      INTEGER nl(3)
      DOUBLE PRECISION amp(3,maxlin),fas(3,maxlin),fre(3,maxlin)
      COMMON /especlineas/amp,fre,fas,nl
      SAVE /especlineas/
* Orbit
      INTEGER np
      DOUBLE PRECISION t(0:nn-1),x(3,0:nn-1),v(3,0:nn-1)
*      COMMON /orbita/t,x,v,np
*      SAVE /orbita/

      LOGICAL prim
      INTEGER i,k
      DOUBLE PRECISION input(2,nn),output(3,maxlin),dpi,tt
      DATA prim/.TRUE./
      SAVE prim,dpi

      IF(prim)THEN
         prim=.FALSE.
         dpi=2d0*ACOS(-1d0)
      ENDIF
      IF(jpan.NE.0)WRITE(*,*)'Computing line spectra...'
* Translate to fmft variables
      tt=t(np-1)+t(1)
      DO k=1,jdim
         DO i=1,np
            input(1,i)=x(k,i-1)
            IF(metodo.EQ.1)THEN
               input(2,i)=v(k,i-1)
            ELSE
               input(2,i)=0d0
            ENDIF
         ENDDO
         CALL fmft(np,input,tt,output)
* Get back to taxon variables
         DO i=1,maxlin
            fre(k,i)=output(1,i)
            amp(k,i)=output(2,i)
            fas(k,i)=MOD(output(3,i)+dpi,dpi)
         ENDDO
      ENDDO
      END

************************************************************************
      SUBROUTINE fmft(ndata,input,tu,output)
      INCLUDE 'taxon.def'
      INTEGER ndata
      DOUBLE PRECISION output(3,maxlin),tu,input(2,ndata)
* Frequency Modified Fourier Transform (Sidlichovsky and Nesvorny 1997,
* Cel. Mech. 65, 137). 
* The vectors input(1,j) and input(2,j), j = 1 ... ndata (ndata must
* be a power of 2), are the input data x(j-1) and v(j-1).
* In the output array, output(1,i), output(2,i) and output(3,i) are
* the i-th frequency, amplitude and phase; tu is the time span.
************************************************************************
      DOUBLE PRECISION fmft_tol,fmft_near
      PARAMETER (fmft_tol=1d-10) ! mft nominal precision
      PARAMETER (fmft_near=0d0)  ! mft overlap exclusion parameter

      LOGICAL prim
      INTEGER nearfreqflag,l,m,j,i,k
      DOUBLE PRECISION pi,xdata(nn),ydata(nn),x(nn),
     & y(nn),powsd(nn),freq(3,maxlin),ampli(3,maxlin),
     & phase(3,maxlin),f(maxlin),a(maxlin),psi(maxlin),
     & q(maxlin,maxlin),alpha(maxlin,maxlin),b(maxlin),centerf,bracket,
     & leftf,rightf,golden,phisqr,fac,xsum,ysum,minfreq,maxfreq
      EXTERNAL phisqr
      DATA prim/.TRUE./
      SAVE prim,pi,minfreq,maxfreq

      DOUBLE PRECISION fftol,twopi
      COMMON /param2/fftol,twopi
      SAVE /param2/

      IF(prim)THEN
         prim=.FALSE.
         pi=ACOS(-1d0)
         twopi=2d0*pi
         minfreq=-pi
         maxfreq=pi
         fftol=fmft_tol
      ENDIF

* Null input
      j=1
      DOWHILE(input(1,j).EQ.0d0.AND.input(2,j).EQ.0d0.AND.j.LT.ndata)
         j=j+1
      ENDDO
      IF(j.EQ.ndata)THEN ! all coordinates zero
         DO k=1,3
            DO i=1,maxlin
               output(k,i)=0d0
            ENDDO
         ENDDO
         RETURN
      ENDIF
* Loops for fmft
      DO l=1,2
         IF(l.EQ.1)THEN
            DO k=1,maxlin
               freq(l,k)=0d0
               ampli(l,k)=0d0
               phase(l,k)=0d0
            ENDDO
* Separate real and imaginary parts
            DO j=1,ndata
               xdata(j)=input(1,j)
               ydata(j)=input(2,j)
            ENDDO
         ELSE
* Generate the quasiperiodic function computed by mft
            DO i=1,ndata
               xdata(i)=0d0
               ydata(i)=0d0
               DO k=1,maxlin
                  xdata(i)=xdata(i)+ampli(l-1,k)*
     &             COS(freq(l-1,k)*(i-1)+phase(l-1,k))
                  ydata(i)=ydata(i)+ampli(l-1,k)*
     &             SIN(freq(l-1,k)*(i-1)+phase(l-1,k))
               ENDDO
            ENDDO
         ENDIF
* Multiply the signal by a window function, store result in x and y
         CALL window(x,y,xdata,ydata,ndata)
* Compute power spectral density using fast Fourier transform
         CALL power(powsd,x,y,ndata)
         IF(l.EQ.1)THEN
* Check if the frequency is in the required range
            centerf=bracket(powsd,ndata)
            DOWHILE(centerf.LT.minfreq.OR.centerf.GT.maxfreq)
* If no, substract it from the signal
               leftf=centerf-twopi/DBLE(ndata)
               rightf=centerf+twopi/DBLE(ndata)
               f(1)=golden(phisqr,leftf,centerf,rightf,x,y,ndata)
               CALL amph(a(1),psi(1),f(1),x,y,ndata)
               DO j=1,ndata
                  xdata(j)=xdata(j)-a(1)*COS(f(1)*(j-1)+psi(1))
                  ydata(j)=xdata(j)-a(1)*SIN(f(1)*(j-1)+psi(1))
               ENDDO
               CALL window(x,y,xdata,ydata,ndata)
               CALL power(powsd,x,y,ndata)
               centerf=bracket(powsd,ndata)
            ENDDO
         ELSE
            centerf=freq(1,1)
         ENDIF
         leftf=centerf-twopi/DBLE(ndata)
         rightf=centerf+twopi/DBLE(ndata)
* Determine the first frequency
         f(1)=golden(phisqr,leftf,centerf,rightf,x,y,ndata)
* Compute amplitude and phase
         CALL amph(a(1),psi(1),f(1),x,y,ndata)
* Substract the first harmonic from the signal
         DO j=1,ndata
            xdata(j)=xdata(j)-a(1)*COS(f(1)*(j-1)+psi(1))
            ydata(j)=ydata(j)-a(1)*SIN(f(1)*(j-1)+psi(1))
         ENDDO
* Here starts the main loop
         q(1,1)=1d0
         alpha(1,1)=1d0 
         DO m=2,maxlin
* Multiply signal by window function
            CALL window(x,y,xdata,ydata,ndata)
* Compute power spectral density using fast Fourier transform
            CALL power(powsd,x,y,ndata)
            IF(l.EQ.1)THEN
               centerf=bracket(powsd,ndata)
               leftf=centerf-twopi/DBLE(ndata)
               rightf=centerf+twopi/DBLE(ndata)
               f(m)=golden(phisqr,leftf,centerf,rightf,x,y,ndata)
* Check whether the new frequency is not too close to any previously
* determined one
               nearfreqflag=0
               DO k=1,m-1
                  IF(ABS(f(m)-f(k)).LT.fmft_near*twopi/DBLE(ndata))
     &             nearfreqflag=1
               ENDDO
* Check if the frequency is in the required range
               DOWHILE(f(m).LT.minfreq.OR.f(m).GT.maxfreq.OR.
     &          nearfreqflag.EQ.1)
* If no, substract it from the signal
                  leftf=centerf-twopi/DBLE(ndata)
                  rightf=centerf+twopi/DBLE(ndata)
                  f(m)=golden(phisqr,leftf,centerf,rightf,x,y,ndata)
                  CALL amph(a(m),psi(m),f(m),x,y,ndata)
                  DO j=1,ndata
                     xdata(j)=xdata(j)-a(m)*COS(f(m)*(j-1)+psi(m))
                     ydata(j)=ydata(j)-a(m)*SIN(f(m)*(j-1)+psi(m))
                  ENDDO
* And recompute the new one
                  CALL window(x,y,xdata,ydata,ndata)
                  CALL power(powsd,x,y,ndata)
                  centerf=bracket(powsd,ndata)
                  leftf=centerf-twopi/DBLE(ndata)
                  rightf=centerf+twopi/DBLE(ndata)
                  f(m)=golden(phisqr,leftf,centerf,rightf,x,y,ndata)
                  nearfreqflag=0
                  DO k=1,m-1
                     IF(ABS(f(m)-f(k)).LT.fmft_near*twopi/DBLE(ndata))
     &                nearfreqflag=1
                  ENDDO
               ENDDO
            ELSE
               centerf=freq(1,m)
               leftf=centerf-twopi/DBLE(ndata)
               rightf=centerf+twopi/DBLE(ndata)
* Determine the next frequency
               f(m)=golden(phisqr,leftf,centerf,rightf,x,y,ndata)
            ENDIF
* Compute its amplitude and phase
            CALL amph(a(m),psi(m),f(m),x,y,ndata)
* Equation (3) in Sidlichovsky and Nesvorny (1997)
            q(m,m)=1d0
            DO j=1,m-1
               fac=(f(m)-f(j))*DBLE(ndata-1)/2d0
               q(m,j)=SIN(fac)/fac*pi**2/(pi**2-fac**2)
               q(j,m)=q(m,j)
            ENDDO
* Equation (17)
            DO k=1,m-1
               b(k)=0d0
               DO j=1,k
                  b(k)=b(k)-alpha(k,j)*q(m,j)
               ENDDO
            ENDDO
* Equation (18)
            alpha(m,m)=1d0
            DO j=1,m-1
               alpha(m,m)=alpha(m,m)-b(j)**2
            ENDDO
            alpha(m,m)=1d0/SQRT(alpha(m,m))
* Equation (19)
            DO k=1,m-1
               alpha(m,k)=0d0
               DO j=k,m-1
                  alpha(m,k)=alpha(m,k)+b(j)*alpha(j,k)
               ENDDO
               alpha(m,k)=alpha(m,m)*alpha(m,k)
            ENDDO
* Equation (22)
            DO i=1,ndata
               xsum=0d0
               ysum=0d0
               DO j=1,m
                  fac=f(j)*(i-1)+(f(m)-f(j))*DBLE(ndata-1)/2d0+psi(m)
                  xsum=xsum+alpha(m,j)*COS(fac)
                  ysum=ysum+alpha(m,j)*SIN(fac)
               ENDDO
               xdata(i)=xdata(i)-alpha(m,m)*a(m)*xsum
               ydata(i)=ydata(i)-alpha(m,m)*a(m)*ysum
            ENDDO
         ENDDO
* Equation (26)
         DO k=1,maxlin
            xsum=0d0
            ysum=0d0
            DO j=k,maxlin
               fac=(f(j)-f(k))*DBLE(ndata-1)/2d0+psi(j)
               xsum=xsum+alpha(j,j)*alpha(j,k)*a(j)*COS(fac)
               ysum=ysum+alpha(j,j)*alpha(j,k)*a(j)*SIN(fac)
            ENDDO
            a(k)=SQRT(xsum**2+ysum**2)
            psi(k)=ATAN2(ysum,xsum)
         ENDDO
* Remember the computed values for the fmft
         DO k=1,maxlin
            freq(l,k)=f(k)
            ampli(l,k)=a(k)
            phase(l,k)=psi(k)
         ENDDO
      ENDDO
* Return the final frequencies, amplitudes and phases
      DO k=1,maxlin
         output(1,k)=freq(1,k)+(freq(1,k)-freq(2,k))
         output(2,k)=ampli(1,k)+(ampli(1,k)-ampli(2,k))
         output(3,k)=phase(1,k)+(phase(1,k)-phase(2,k))
      ENDDO
* Sort the frequencies in decreasing order of amplitude
      CALL dsort(output)
      DO k=1,maxlin
         output(1,k)=output(1,k)/(twopi*tu/DBLE(ndata))
         output(3,k)=MOD(output(3,k)+twopi,twopi)
      ENDDO
      END  

************************************************************************
      SUBROUTINE window(x,y,xdata,ydata,ndata)
      INTEGER ndata
      DOUBLE PRECISION xdata(ndata),ydata(ndata),x(ndata),y(ndata)
* Multiplies data by a window function
************************************************************************
      INTEGER j
      DOUBLE PRECISION win

      DOUBLE PRECISION fftol,twopi
      COMMON /param2/fftol,twopi
      SAVE /param2/

      DO j=1,ndata
         win=twopi*DBLE(j-1)/DBLE(ndata-1)
         win=(1d0-COS(win))/2d0
         x(j)=xdata(j)*win
         y(j)=ydata(j)*win
      ENDDO
      END

************************************************************************
      SUBROUTINE power(powsd,x,y,ndata)
      INTEGER ndata
      DOUBLE PRECISION x(ndata),y(ndata),powsd(ndata)
* Rearranges data for the fast Fourier transform, calls fft and returns
* power spectral density.
************************************************************************
      INCLUDE 'taxon.def'

      INTEGER j
      DOUBLE PRECISION z(2*nn)

      DO j=1,ndata
         z(2*j-1)=x(j)
         z(2*j)=y(j)
      ENDDO
      CALL four11(z,ndata,1)
      DO j=1,ndata
         powsd(j)=z(2*j-1)**2+z(2*j)**2
      ENDDO
      END
 
************************************************************************
      SUBROUTINE four11(dato,nn,isign)
      INTEGER nn,isign
      DOUBLE PRECISION dato(2*nn)
* dato[1..2*nn] replaces by DFS, nn must be a power of 2
************************************************************************
      INTEGER n,mmax,m,j,istep,i
      DOUBLE PRECISION wtemp,wr,wpr,wpi,wi,theta,tempr,tempi,aux

      DOUBLE PRECISION fftol,twopi
      COMMON /param2/fftol,twopi
      SAVE /param2/

      n=nn*2
      j=1
* Bit-reversal section
      DO i=1,n-1,2
         IF(j.GT.i)THEN
            aux=dato(j)
            dato(j)=dato(i)
            dato(i)=aux
            aux=dato(j+1)
            dato(j+1)=dato(i+1)
            dato(i+1)=aux
         ENDIF
         m=n/2
         DOWHILE(m.GE.2.AND.j.GT.m)
            j=j-m
            m=m/2
         ENDDO
         j=j+m
      ENDDO
* Danielson-Lanczos section
      mmax=2
      DOWHILE(n.GT.mmax) ! outer ln nn loop
         istep=mmax*2
         theta=isign*twopi/DBLE(mmax) ! initialize
         wtemp=SIN(0.5d0*theta)
         wpr=-2d0*wtemp**2
         wpi=SIN(theta)
         wr=1d0
         wi=0d0
         DO m=1,mmax-1,2 ! two inner loops
            DO i=m,n,istep
               j=i+mmax ! D-L formula
               tempr=wr*dato(j)-wi*dato(j+1)
               tempi=wr*dato(j+1)+wi*dato(j)
               dato(j)=dato(i)-tempr
               dato(j+1)=dato(i+1)-tempi
               dato(i)=dato(i)+tempr
               dato(i+1)=dato(i+1)+tempi
            ENDDO
            wtemp=wr
            wr=wtemp*wpr-wi*wpi+wr ! trig. recurrence
            wi=wi*wpr+wtemp*wpi+wi
         ENDDO
         mmax=istep
      ENDDO
      END

************************************************************************
      FUNCTION bracket(powsd,ndata)
      INTEGER ndata
      DOUBLE PRECISION bracket,powsd(ndata)
* Finds the maximum of the power spectral density
************************************************************************
      INTEGER j,maxj
      DOUBLE PRECISION freq,maxpow

      DOUBLE PRECISION fftol,twopi
      COMMON /param2/fftol,twopi
      SAVE /param2/

      maxj=0
      maxpow=0d0
      DO j=2,ndata/2-2
         IF(powsd(j).GT.powsd(j-1).AND.powsd(j).GT.powsd(j+1))THEN
            IF(powsd(j).GT.maxpow)THEN
               maxj=j
               maxpow=powsd(j)
            ENDIF
         ENDIF
      ENDDO
      DO j=ndata/2+2,ndata-1
         IF(powsd(j).GT.powsd(j-1).AND.powsd(j).GT.powsd(j+1))THEN
            IF(powsd(j).GT.maxpow)THEN
               maxj=j
               maxpow=powsd(j)
            ENDIF
         ENDIF
      ENDDO
      IF(powsd(1).GT.powsd(2).AND.powsd(1).GT.powsd(ndata))THEN
         IF(powsd(1).GT.maxpow)THEN
            maxj=j
            maxpow=powsd(j)
         ENDIF
      ENDIF
* This should not happen because we eliminated null coordinates
*     IF(maxpow.EQ.0d0)STOP 'bracket: DFT has no maximum'
* Negative signs and twopi compensate for the Numerical Recipes 
* definition of the DFT
      IF(maxj.LT.ndata/2)freq=-(maxj-1)
      IF(maxj.GT.ndata/2)freq=-(maxj-ndata-1)
      bracket=twopi*freq/DBLE(ndata)
      END

************************************************************************
      FUNCTION golden(f,ax,bx,cx,xdata,ydata,n)
      INTEGER n
      DOUBLE PRECISION golden,f,ax,bx,cx,xdata(n),ydata(n)
      EXTERNAL f
* Calculates the maximum of a function bracketed by ax, bx and cx
************************************************************************
      DOUBLE PRECISION gold_r,gold_c
      PARAMETER (gold_r=0.61803399d0)
      PARAMETER (gold_c=1d0-gold_r)

      DOUBLE PRECISION fftol,twopi
      COMMON /param2/fftol,twopi
      SAVE /param2/

      DOUBLE PRECISION f1,f2,x0,x1,x2,x3

      x0=ax
      x3=cx
      IF(ABS(cx-bx).GT.ABS(bx-ax))THEN
         x1=bx
         x2=bx+gold_c*(cx-bx)
      ELSE
         x2=bx
         x1=bx-gold_c*(bx-ax)
      ENDIF
      f1=f(x1,xdata,ydata,n)
      f2=f(x2,xdata,ydata,n)
c      DOWHILE(ABS(x3-x0).GT.fftol*(ABS(x1)+ABS(x2)))
      DOWHILE(ABS(x3-x0).GT.fftol*(ABS(x1)+ABS(x2)).AND.
     & ABS(x3-x0).GT.fftol)
         IF(f2.GT.f1)THEN
            x0=x1
            x1=x2
            x2=gold_r*x1+gold_c*x3
            f1=f2
            f2=f(x2,xdata,ydata,n)
         ELSE
            x3=x2
            x2=x1
            x1=gold_r*x2+gold_c*x0
            f2=f1
            f1=f(x1,xdata,ydata,n)
         ENDIF
      ENDDO
      IF(f1.GT.f2)THEN
         golden=x1
      ELSE
         golden=x2
      ENDIF
      END

************************************************************************
      SUBROUTINE amph(amp,phase,freq,xdata,ydata,ndata)
      INTEGER ndata
      DOUBLE PRECISION amp,phase,freq,xdata(ndata),ydata(ndata)
* Calculates the amplitude and phase
************************************************************************
      DOUBLE PRECISION xphi,yphi

      xphi=0d0
      yphi=0d0
      CALL phifun(xphi,yphi,freq,xdata,ydata,ndata)
      amp=SQRT(xphi**2+yphi**2)
      phase=ATAN2(yphi,xphi)
      END
  
************************************************************************
      FUNCTION phisqr(freq,xdata,ydata,ndata)
      INTEGER ndata
      DOUBLE PRECISION phisqr,freq,xdata(ndata),ydata(ndata)
* Computes a square power of the function phi
************************************************************************
      INCLUDE 'taxon.def'

      INTEGER i,ndatal
      DOUBLE PRECISION xphi,yphi,freql,xdatal(nn),ydatal(nn)

* To local variables
      ndatal=ndata
      freql=freq
      DO i=1,ndata
         xdatal(i)=xdata(i)
         ydatal(i)=ydata(i)
      ENDDO
      xphi=0d0
      yphi=0d0
      CALL phifun(xphi,yphi,freql,xdatal,ydatal,ndatal)
      phisqr=xphi**2+yphi**2
      END

************************************************************************
      SUBROUTINE phifun(xphi,yphi,freq,xdata,ydata,n)
      INTEGER n
      DOUBLE PRECISION xphi,yphi,freq,xdata(n),ydata(n)
* Computes the function phi 
************************************************************************
      INCLUDE 'taxon.def'

      INTEGER i,j,na
      DOUBLE PRECISION c,s,xdata2(nn),ydata2(nn)

      xdata2(1)=xdata(1)/2d0
      ydata2(1)=ydata(1)/2d0
      xdata2(n)=xdata(n)/2d0
      ydata2(n)=ydata(n)/2d0
      DO i=2,n-1
         xdata2(i)=xdata(i)
         ydata2(i)=ydata(i)
      ENDDO
      na=n
      DOWHILE(na.NE.1)
         na=na/2
         c=COS(-na*freq)
         s=SIN(-na*freq)
         DO i=1,na
            j=i+na
            xdata2(i)=xdata2(i)+c*xdata2(j)-s*ydata2(j)
            ydata2(i)=ydata2(i)+c*ydata2(j)+s*xdata2(j)
         ENDDO
      ENDDO
      xphi=2d0*xdata2(1)/DBLE(n-1)
      yphi=2d0*ydata2(1)/DBLE(n-1)
      END

************************************************************************
      SUBROUTINE dsort(output)
      INCLUDE 'taxon.def'
      DOUBLE PRECISION output(3,maxlin)
* Sorting procedure
************************************************************************
      INTEGER j,iwksp(maxlin),n2
      DOUBLE PRECISION wksp(maxlin)

      n2=maxlin+1

      DO j=1,maxlin
         wksp(j)=output(2,j)
      ENDDO
      CALL quicki8(maxlin,wksp,iwksp)
      DO j=1,maxlin
         output(2,j)=wksp(n2-j)
      ENDDO
      DO j=1,maxlin
         wksp(j)=output(1,j)
      ENDDO
      DO j=1,maxlin
         output(1,j)=wksp(iwksp(n2-j))
      ENDDO
      DO j=1,maxlin
         wksp(j)=output(3,j)
      ENDDO
      DO j=1,maxlin
         output(3,j)=wksp(iwksp(n2-j))
      ENDDO
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                       MISCELANEOUS SECTION                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

************************************************************************
      FUNCTION is0(x1,x2)
      LOGICAL is0
      DOUBLE PRECISION x1,x2
* It determines whether x1+x2 is zero, by means of:
*  |x1+x2| / SQRT(x1**2+x2**2)) < eps
* eps: parameter to set the precision
************************************************************************
      DOUBLE PRECISION eps
      COMMON /epsilon/eps
      SAVE /epsilon/

      is0=.FALSE.
* Comparing with zero
      IF((x1.EQ.0d0.AND.ABS(x2).LT.eps).OR.
     & (x2.EQ.0d0.AND.ABS(x1).LT.eps))THEN
         is0=.TRUE.
         RETURN
      ENDIF
* Any other comparison
      IF(ABS(x1+x2)/SQRT(x1**2+x2**2).LT.eps)is0=.TRUE.
      END

************************************************************************
      SUBROUTINE quicki8(n,x,ind)
      INTEGER n,ind(n)
      DOUBLE PRECISION x(n)
* It sorts an array (smallest to greatest) and another one of indexes
* according to it. 
* Input:  n   - number of elements to be sorted
*         x   - array to be sorted
*         ind - array of indexes to be sorted according to x
* Output: x   - sorted vector
*         ind - ind(k) is the former index of x(k); another array y 
*               will have y(ind(k)) corresp. to x(k)
* If n is less than the real dimension of x, it sorts the first n 
* elements, leaving the rest untouched.
************************************************************************
      LOGICAL flag,fleg,flig
      INTEGER auxx,auxx2,i,is,ita(2,60),j,l,m,r
      DOUBLE PRECISION aux,aux2

      m=9
* Initialization of the indirect addressing
      DO i=1,n
         ind(i)=i
      ENDDO
* Quicksort algorithm (Knuth)
*---Q1------------------------------------------------------------------
      l=1 
      r=n 
      is=0
*---Q2------------------------------------------------------------------
      fleg=.TRUE.
      DOWHILE(fleg)
         fleg=.FALSE.
         i=l 
         j=r+1 
         aux=x(l)
         auxx=ind(l)
*---Q3------------------------------------------------------------------
         flag=.TRUE.
         DOWHILE(flag)
            i=i+1  
            flig=.TRUE.
            DOWHILE(i.NE.n+1.AND.flig)
               IF(x(i).LT.aux)THEN
                  i=i+1 
               ELSE
                  flig=.FALSE.
               ENDIF
            ENDDO
*---Q4------------------------------------------------------------------
            j=j-1 
            DOWHILE(x(j).GT.aux)
               j=j-1
            ENDDO
*---Q6------------------------------------------------------------------
            IF(j.GT.i)THEN 
               aux2=x(i)  
               x(i)=x(j) 
               x(j)=aux2
               auxx2=ind(i)
               ind(i)=ind(j)
               ind(j)=auxx2
            ELSE
               flag=.FALSE.
            ENDIF
         ENDDO
*---Q5------------------------------------------------------------------
         x(l)=x(j) 
         x(j)=aux
         ind(l)=ind(j)
         ind(j)=auxx
*---Q7------------------------------------------------------------------
         IF(r-j.GE.j-l)THEN
            IF(j-l.GT.m)THEN 
               is=is+1 
               ita(1,is)=j+1 
               ita(2,is)=r
               r=j-1 
               fleg=.TRUE.
            ELSE 
               IF(r-j.GT.m)THEN 
                  l=j+1 
                  fleg=.TRUE.
*---Q8------------------------------------------------------------------
               ELSEIF(is.NE.0)THEN
                  l=ita(1,is) 
                  r=ita(2,is) 
                  is=is-1 
                  fleg=.TRUE.
*-----------------------------------------------------------------------
               ENDIF
            ENDIF
         ELSE
            IF(r-j.GT.m)THEN 
               is=is+1 
               ita(1,is)=l 
               ita(2,is)=j-1
               l=j+1 
               fleg=.TRUE.
            ELSE 
               IF(j-l.GT.m)THEN 
                  r=j-1 
                  fleg=.TRUE.
*---Q8------------------------------------------------------------------
               ELSEIF(is.NE.0)THEN
                  l=ita(1,is) 
                  r=ita(2,is) 
                  is=is-1 
                  fleg=.TRUE.
*-----------------------------------------------------------------------
               ENDIF
            ENDIF
         ENDIF
      ENDDO
*---Q9------------------------------------------------------------------
      DO j=2,n
         IF(x(j-1).GT.x(j))THEN
            i=j-1
            aux=x(j)
            auxx=ind(j)
            x(i+1)=x(i)
            ind(i+1)=ind(i)
            i=i-1
            fleg=.TRUE.
            DOWHILE(i.NE.0.AND.fleg)
               IF(aux.LT.x(i))THEN
                  x(i+1)=x(i)
                  ind(i+1)=ind(i)
                  i=i-1 
               ELSE
                  fleg=.FALSE.
               ENDIF
            ENDDO
            x(i+1)=aux
            ind(i+1)=auxx
         ENDIF
      ENDDO
*---End quicksort-------------------------------------------------------
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         INPUT/OUTPUT SECTION                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

************************************************************************
      SUBROUTINE  entrada(ju,np,t,x,v,jsub,jdim,jcla,jcl,jpan,jlin,
     \                   jcom,arch)
      INTEGER ju
* Input of data
************************************************************************
      INCLUDE 'taxon.def'
* Orbit
      INTEGER np
      DOUBLE PRECISION t(0:nn-1),x(3,0:nn-1),v(3,0:nn-1)
*      COMMON /orbita/t,x,v,np
*      SAVE /orbita/
* Input
      INTEGER jsub,jdim,jcla,jcl,jpan,jlin,jcom
      CHARACTER arch*30
*      COMMON /input/jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch
      
* Equality of frequencies
      DOUBLE PRECISION eps
      COMMON /epsilon/eps
      SAVE /epsilon/
* Line spectra
      INTEGER nl(3)
      DOUBLE PRECISION amp(3,maxlin),fas(3,maxlin),fre(3,maxlin)
      COMMON /especlineas/amp,fre,fas,nl
      SAVE /especlineas/

      LOGICAL dos
      INTEGER i,k
      DOUBLE PRECISION ttot

* Verify if taxon was called as a subroutine
      IF(jsub.EQ.0)THEN  ! taxon is a program; read from file
         OPEN(50,FILE='taxon.par')
         READ(50,'(a)')arch
         READ(50,*)jdim,jcla,jcl,jpan,jlin
         CLOSE(50)
         jcom=0
      ENDIF
      IF(jdim.NE.2.AND.jdim.NE.3)STOP'ERROR: jdim must be 2 or 3'
* Non-blank characters of arch, except the extension, if any.
      ju=INDEX(arch,' ')-1
      IF(INDEX(arch,'.').NE.0)ju=INDEX(arch,'.')-1
* Input orbit
c      IF(jcom.EQ.0)THEN ! Read from file
c         IF(jpan.NE.0)WRITE(*,*)'Reading input data...'
c         CALL orbitin(jdim,arch)
c      ENDIF
* Check np = 2**n
      dos=.FALSE.
      DO i=1,30
         IF(np.EQ.2**i)dos=.TRUE.
      ENDDO
      IF(.NOT.dos)STOP'ERROR: np must be a power of 2.'
      IF(np.GT.nn)STOP'ERROR: np > nn.'
* In case of a 2D orbit with no z-input
      IF(jdim.EQ.2)THEN
         DO i=0,nn-1
            x(3,i)=0d0
         ENDDO
      ENDIF
* Time of integration: N*Delta (one more Delta than t(np-1))
      ttot=(t(np-1)-t(0))+(t(1)-t(0))
* eps as a function of the integration period (1/T = fj(1,1) = width of
* a Fourier slot)
      IF(jdim.EQ.2)THEN
         eps=epsf2/ttot
      ELSE
         eps=epsf3/ttot
      ENDIF
      eps=MIN(eps,epsmin)
c      WRITE(*,*)'epsmin = ',epsmin
c      WRITE(*,*)'eps = ',eps
* Initialize the spectra
      DO i=1,maxlin
         DO k=1,3
             fre(k,i)=0d0
             amp(k,i)=0d0
             fas(k,i)=0d0
         ENDDO
      ENDDO
      END

************************************************************************
      SUBROUTINE salida(ju,np,t,x,v,jsub,jdim,jcla,jcl,jpan,jlin,
     \                  jcom,arch,code,peri)
      INTEGER ju
* Output to file.
************************************************************************
      INCLUDE 'taxon.def'
* Input
      INTEGER jsub,jdim,jcla,jcl,jpan,jlin,jcom
      CHARACTER arch*30
*      COMMON /input/jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch
* Output
      INTEGER code,lic1(3),lic2(3),licf(3),peri
      DOUBLE PRECISION ffund(4)
      COMMON /results/ffund,lic1,lic2,licf
      SAVE /results/
* Orbit
      INTEGER np
      DOUBLE PRECISION t(0:nn-1),x(3,0:nn-1),v(3,0:nn-1)
*      COMMON /orbita/t,x,v,np
*      SAVE /orbita/
* Line spectra
      INTEGER nl(3)
      DOUBLE PRECISION amp(3,maxlin),fas(3,maxlin),fre(3,maxlin)
      COMMON /especlineas/amp,fre,fas,nl
      SAVE /especlineas/

      LOGICAL prim
      INTEGER j,k,l,nff,nres,u
      DOUBLE PRECISION dospi,pi,ttot
      CHARACTER pref(3)*1,tipo*15
      DATA prim/.TRUE./
      SAVE prim,pi,dospi,pref

      IF(prim)THEN
         prim=.FALSE.
         pi=ACOS(-1d0)
         dospi=2d0*pi
         pref(1)='x'
         pref(2)='y'
         pref(3)='z'
      ENDIF
* Output of the raw spectra. Each line consists in three records: 
*          frequency, 0, 0
*          frequency, amplitude, phase in degrees
*          [blank line]
* Thus, these spectra can be plotted using the convention that a blank
* record signals a new line to be plotted.
      IF(jlin.NE.0)THEN
         DO k=1,jdim
            OPEN(61,FILE=arch(1:ju)//'.li'//pref(k))
            DO l=1,nl(k)
               WRITE(61,*)REAL(fre(k,l)),0.,0.
               WRITE(61,*)REAL(fre(k,l)),REAL(amp(k,l)),
     &            REAL(180d0*fas(k,l)/pi)
               WRITE(61,*)
            ENDDO
            CLOSE(61)
         ENDDO
* Output of the spectra in logarithmic scale. Each line consists in
* three records: frequency, -10
*                frequency, LOG(amplitude)
*                [blank line]
* Thus, these spectra can be plotted using the convention that a blank
* record signals a new line to be plotted. There is also a third
* value in the first two records: the frequency in radians per unit
* time, in case these values are wanted. 
         DO k=1,jdim
            OPEN(62,FILE=arch(1:ju)//'.ll'//pref(k))
            DO l=1,nl(k)
               IF(amp(k,l).GT.0d0)THEN
                  WRITE(62,*)REAL(fre(k,l)),-10.,REAL(dospi*fre(k,l))
                  WRITE(62,*)REAL(fre(k,l)),REAL(LOG10(amp(k,l))),
     &               REAL(dospi*fre(k,l))
                  WRITE(62,*)
               ENDIF
            ENDDO
            CLOSE(62)
         ENDDO
      ENDIF
* Orbital periods
      ttot=(t(np-1)-t(0))+(t(1)-t(0))
      IF(nl(1).NE.0)THEN
         peri=NINT(ABS(fre(1,1))*ttot)
      ELSEIF(nl(2).NE.0)THEN
         peri=NINT(ABS(fre(2,1))*ttot)
      ELSE
         peri=NINT(ABS(fre(3,1))*ttot)
      ENDIF
* Number of base frequencies, of resonances, and morphology of the orbit
      nff=code/100
      nres=(code-100*nff)/10
      u=MOD(code,10)
* Name of the orbit
      tipo='               '
      IF(nff.EQ.4)THEN
         tipo='irregular'
      ELSEIF(nff.EQ.3)THEN
         IF(jdim.EQ.3)THEN
            tipo='open'
         ELSE
            tipo='irregular'
         ENDIF
      ELSEIF(nff.EQ.2)THEN
         IF(jdim.EQ.3)THEN
            tipo='thin'
         ELSE
            tipo='open'
         ENDIF
      ELSEIF(nff.EQ.1)THEN
         tipo='closed'
      ELSE
         tipo='not classified'
      ENDIF
      IF(tipo(1:3).NE.'irr'.AND.tipo(1:3).NE.'not')THEN
         IF(jdim.EQ.2)THEN
            IF(u.EQ.0)THEN
               tipo(8:11)='box'
            ELSE
               tipo(8:12)='loop'
            ENDIF
         ELSE
            IF(u.EQ.0)THEN
               tipo(8:15)='box     '
            ELSEIF(u.EQ.1)THEN
               tipo(8:15)='x tube  '
               CALL tubox(tipo(15:15),np,t,x,v)
            ELSEIF(u.EQ.2)THEN
               tipo(8:15)='y tube  '
            ELSE
               tipo(8:15)='z tube  '
            ENDIF
         ENDIF
      ENDIF
* Output to screen
      IF(jpan.NE.0)THEN
         WRITE(*,'(a)')'--------------'
         WRITE(*,'(a)')'Classification'
         WRITE(*,'(a)')'--------------'
         WRITE(*,*)'Code = ',code
         WRITE(*,*)'Name = ',tipo
         WRITE(*,*)nff,' <-- # of fundamental frequencies'
         DO j=1,nff
            IF(j.EQ.1)THEN
               WRITE(*,*)ffund(1),' <-- fundamental freq. [cycles/u.t.]'
            ELSE
               WRITE(*,*)ffund(j)
            ENDIF
         ENDDO
         WRITE(*,*)nres,' <-- # of resonances'
         IF(nres.GE.1)WRITE(*,*)'( ',lic1(1),':',lic1(2),':',lic1(3),')'
         IF(nres.EQ.2)WRITE(*,*)'( ',lic2(1),':',lic2(2),':',lic2(3),')'
         WRITE(*,*)peri,' <-- Orbital periods'
      ENDIF
* Output to file
      IF(jcla.NE.0)THEN
         OPEN(64,FILE=arch(1:ju)//'.cla')
         WRITE(64,*)'Classification code = ',code
         WRITE(64,*)'Name = ',tipo
         WRITE(64,*)'# of fundamental frequencies = ',nff
         DO j=1,nff
            IF(j.EQ.1)THEN
               WRITE(64,*)'Fundamental freq. [cycles/u.t.] = ',ffund(1)
            ELSE
               WRITE(64,*)'                                 ',ffund(j)
            ENDIF
         ENDDO
         WRITE(64,*)'# of resonances = ',nres
         IF(nres.GE.1)WRITE(64,*)'( ',lic1(1),':',lic1(2),':',lic1(3),
     &      ')'
         IF(nres.EQ.2)WRITE(64,*)'( ',lic2(1),':',lic2(2),':',lic2(3),
     &      ')'
         WRITE(64,*)'Orbital periods = ',peri
         WRITE(64,*)
         CLOSE(64)
      ENDIF
* Simplified output, with only relevant results and without banners.
* Ideal for, e.g., massive orbit classification. 
      IF(jcl.NE.0)THEN
         OPEN(65,FILE=arch(1:ju)//'.cl')
         DOWHILE(.TRUE.)
            READ(65,'(a)',END=1)
         ENDDO
1        BACKSPACE(65)
         IF(jdim.EQ.3)THEN
            WRITE(65,'(8i6)')code,(lic1(k),k=1,3),(lic2(k),k=1,3)
         ELSE
            WRITE(65,'(3i6)')code,(lic1(k),k=1,2)
         ENDIF
         CLOSE(65)
      ENDIF
* Warnings about orbital periods
      IF(jpan.NE.0)THEN
         IF(peri.LE.80)WRITE(*,'(/a/)')'WARNING: less than 80 periods'
         IF(peri.GT.200)WRITE(*,'(/a/)')'WARNING: more than 200 periods'
      ENDIF
      END

************************************************************************
c      SUBROUTINE orbitin(jdim,arch)
c      INTEGER jdim
c      CHARACTER arch*30
* Orbit input.
* Please re-write this routine as desired in order to conform with the
* format of your orbit file. 
*
* Details to take into account:
* a) You may take the variable "arch" as the name of the file. It is the
*    name already read by the input routine from the file "taxon.par".
* b) Times must be stored in the array DOUBLE PRECISION t(0:nn-1). Here 
*    "nn" is a parameter set in "taxon.def"; the actual number of points
*    (np) should be less than or equal to nn.
* c) Positions and velocities, which should be equidistant in time, must
*    be stored in the arrays DOUBLE PRECISION x(3,0:nn-1),v(3,0:nn-1).
*    The number of points must be a power of two.
* d) The routine must set the variable "np" to the actual number of
*    points of the orbit.
************************************************************************
c      INCLUDE 'taxon.def'
* Orbit
c      INTEGER np
c      DOUBLE PRECISION t(0:nn-1),x(3,0:nn-1),v(3,0:nn-1)
c      COMMON /orbita/t,x,v,np
c      SAVE /orbita/

c      INTEGER i

*-----------------------------------------------------------------------
* Write here the input routine------------------------------------------
* This is an example. "arch" is the string read from taxon.par.

c      INTEGER k
* 2D orbit
c      OPEN(51,FILE=arch,STATUS='OLD',ERR=2)
c      DO i=0,nn-1
c         IF(jdim.EQ.2)THEN
c            READ(51,*,ERR=2)t(i),(x(k,i),k=1,2),(v(k,i),k=1,2)
c         ELSE
c            READ(51,*,ERR=2)t(i),(x(k,i),k=1,3),(v(k,i),k=1,3)
c         ENDIF
c      ENDDO
c      CLOSE(51)

c      np=nn

* This is another example, reading from a binary file
*      OPEN(51,FILE=arch,FORM='UNFORMATTED',STATUS='OLD',ERR=2)
*      READ(51,END=1,ERR=2)np,(t(i),i=0,np-1),((x(k,i),i=0,np-1),k=1,3),
*     &   ((v(k,i),i=0,np-1),k=1,3)
*      CLOSE(51)

* End of the input routine----------------------------------------------
*-----------------------------------------------------------------------

c      RETURN
c2     STOP'ERROR: some error while opening/reading the orbit.'
c      END
