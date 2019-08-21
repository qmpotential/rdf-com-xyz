      program rdf_com

C     This program genrerates COM-COM radial distribution using xyz file
C
C     Copyright (C) 2019  Oleg N. Starovoytov
C
C     Last time updated: Wed Aug 21 12:10:04 JST 2019
C
C     This program is free software; you can redistribute it and/or
C     modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2
C     of the License, or (at your option) any later version.
C
C     https://www.gnu.org/licenses/quick-guide-gplv3.html
C  
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.

      implicit none

      integer n
      parameter (n=1000000)

      integer   natm, tatm
      integer   i_unit, j_unit, k_unit
      integer   i_type,  mol_num(n), atm_num(n)
      integer   i_count, mol_count,  atm_count, mol_type
      integer   i, j, k, l, m, o, p, q, s, t, u, v, w
      integer   nfr, num, ig, nbin, g(50,n), i_rdf

      real*4    NA, AN  
      real*4    x(n), y(n), z(n)
      real*4    boxAX,   boxAY,   boxAZ,  volume
      real*4    atm_mass, sys_mass

      real*4    cn(50,n), gr(50,n)
      real*4    massX(50,n), massY(50,n), massZ(50,n), massM(50,n)
      real*4    xx, yy, zz, r, rdb, bin, vbin, rho, nid, nrho(50)

      character*3 ext
      character*4 atm_name
      character*72 file_name
C
      i_unit = 67 
      j_unit = 68 
      k_unit = 69 
C
      ext    = 'RDF'
C
          NA = 6.02214076E+23
          AN = 1E-24

C     Redaing 'rdf.par' parameter file 

      open (i_unit,file='rdf.par',status='old')

      read(i_unit,*)
      read(i_unit,*)  mol_type
      read(i_unit,*)
      read(i_unit,*)  (mol_num(i), i = 1, mol_type)
      read(i_unit,*)
      read(i_unit,*)  (atm_num(j), j = 1, mol_type)
      read(i_unit,*)
      read(i_unit,*)  rdb, nbin
      read(i_unit,*)
      read(i_unit,*)  boxAX, boxAY, boxAZ

      close(i_unit)

C     Initializing number of bins

      do  i = 1, mol_type

           do  j = 1, nbin
          g(i,j) = 0
         cn(i,j) = 0.0d0
         gr(i,j) = 0.0d0
           enddo

      enddo
C
      bin   =  (rdb/nbin)
C
      write(6,*) '-----------------------------------------------------'
      write(6,*) 'Distance, Å   :', rdb
      write(6,*) 'Number of bins:', nbin
      write(6,*) 'Bin size, Å   :', bin
      write(6,*) '-----------------------------------------------------'

C     Calculate total number of atoms

      natm = 0

      do k = 1, mol_type
         natm = natm + int(real(mol_num(k)) * real(atm_num(k)))
      enddo

C     Calculate volume of simulation box

      volume   =    0.0d0
      volume   =    boxAX * boxAY * boxAZ
C
      write(6,*) 'Mol types    :',  mol_type
      write(6,*) 'Num of mols  :', (mol_num(i), i = 1, mol_type)
      write(6,*) 'Num of atoms :', (atm_num(j), j = 1, mol_type)
      write(6,*) 'Total N atoms:',  natm
      write(6,*) '-----------------------------------------------------'

C     Input xyz file

      write(6,*) 'Please type the name of your XYZ file.'
       read(*,*)  file_name
      write(6,*) '-----------------------------------------------------'

      open (j_unit, file=file_name, status='unknown')

C     Initializing number of frames

      nfr  = 0

   34 continue 

C     Start reading xyz file

      read(j_unit,*,end=1289) tatm
      read(j_unit,*,end=1289)

C     Initializing frame count

              nfr = nfr + 1

C     Number of atoms in rdf.par file should be equal to the number in xyz file
     
      if (tatm.ne.natm) write(6,*) 'Total number of atoms is wrong! ';  exit 

C     Initialize parameters for calculatin of radial distributions

         sys_mass = 0.0d0

            i_rdf = 0
           i_type = 1
          i_count = 1
        mol_count = 1
        atm_count = 0

C     Start reading atomic positions

      do  l = 1, tatm

            read(j_unit,*,end=1289) atm_name, x(l), y(l), z(l)

                             atm_count = atm_count + 1

                         if (atm_count.gt.atm_num(i_type)) then 
                               i_count =   i_count + 1
                             mol_count = mol_count + 1
                             atm_count = 1 
                         endif

                         if (mol_count.gt.mol_num(i_type)) then 
                                i_type = i_type  + 1
                             mol_count = 1
                         endif

C     Assigning atomic mass and calculating COM coordinates

       call assign_name(atm_name, atm_mass)

            massX( i_type, mol_count ) = 
     &      massX( i_type, mol_count ) + atm_mass * x(l)

            massY( i_type, mol_count ) =
     &      massY( i_type, mol_count ) + atm_mass * y(l)

            massZ( i_type, mol_count ) =
     &      massZ( i_type, mol_count ) + atm_mass * z(l)

            massM( i_type, mol_count ) = 
     &      massM( i_type, mol_count ) + atm_mass

             sys_mass = sys_mass + atm_mass

      enddo 

C     Finilizing COM coordinate calculations 

                 do m = 1, i_type
                 do o = 1, mol_num(m)

           massX(m,o) = massX(m,o) / massM(m,o)
           massY(m,o) = massY(m,o) / massM(m,o)
           massZ(m,o) = massZ(m,o) / massM(m,o)

                 enddo
                 enddo

C     Start calculations of radial distribution

                 do s = 1, i_type
                 do t = 1, i_type

                 if (t.eq.s.or.t.gt.s) then

                i_rdf = i_rdf + 1

                 do p = 1, mol_num(s)
                 do q = 1, mol_num(t)

                   xx = massX(s,p) - massX(t,q)
                   yy = massY(s,p) - massY(t,q)
                   zz = massZ(s,p) - massZ(t,q)

C     Periodic boundary condition

                   xx = xx - boxAX * nint( xx / boxAX )
                   yy = yy - boxAY * nint( yy / boxAY )
                   zz = zz - boxAZ * nint( zz / boxAZ )

C     COM-COM distance

                   r  = sqrt( xx * xx + yy * yy + zz * zz)

                   if ( r .lt. boxAX / 2.0d0 ) then
                   if ( r .lt. boxAY / 2.0d0 ) then
                   if ( r .lt. boxAZ / 2.0d0 ) then

                  ig  = int( r / bin )
          
          if (ig.ne.0) g(i_rdf,ig) = g(i_rdf,ig) + 1
          
                   endif
                   endif
                   endif
C
                 enddo
                 enddo
                 
                 endif

                 enddo
                 enddo

C     Clear up COM coordinates for the next frame

                 do w = 1, i_type
                 do u = 1, mol_num(w)

                      massX(w,u) = 0.0d0
                      massY(w,u) = 0.0d0
                      massZ(w,u) = 0.0d0
                      massM(w,u) = 0.0d0

                 enddo
                 enddo
C
      goto 34

 1289 continue

C     Calculating the density of a system   

      rho      = (sys_mass / NA) / (volume * AN)

C     Print average results of calculations

      write(6,*) 'Density,     g/cm3:', rho
      write(6,*) 'Average Volume, Å3:', volume
      write(6,*) 'Number of frames  :', nfr
      write(6,*) 'Average box size  :', boxAX, boxAY, boxAZ
      write(6,*) 'System mass       :', sys_mass
      write(6,*) '-----------------------------------------------------'

C     Initializing file for RDF output

      call parse_name(file_name, ext, k_unit)

C     Calculating radial distributions

      call make_rdf(bin, nbin, i_type, i_rdf, nfr, volume, mol_num,
     &                                            g, gr, cn)

C     Print results of radial calculations to a file

      call print_rdf(bin, nbin, i_rdf, gr, cn, k_unit)

      write(6,*) 'CLOSING FILE:      ', file_name
      close(j_unit)
      close(k_unit)
C
      end

C     ------------------------------------------------------------------ 
C                                                 Parsing input filename

      subroutine parse_name(file_name, ext, k_unit)

      integer       i, k_unit, lth
      character*(*) file_name, ext

      lth = LEN(TRIM(file_name))

      do i=1, lth
      if (file_name(i:i) == '.') then
      open(k_unit, file = file_name(1:i)//ext, status='unknown')
      write(6,*) 'OPENING FILE:  ', file_name(1:i)//ext
      endif
      enddo
C
      return
C
      end

C     ------------------------------------------------------------------ 
C                                                               Make rdf

      subroutine make_rdf(bin, nbin, i_type, i_rdf, nfr, volume, 
     &                                               mol_num, g, gr, cn)

C     Frenkel, D. and Berend, S.; Understanding MD simulations: 
C     from algorithms to applications; Chapter 4: Molecular
C     dynamics simulations; Academic Press: Elsevier 2007; pp 85-86.

      integer n
      parameter (n=1000000)

      integer    i, j, k, nbin, i_type, i_rdf, nfr 
      integer    mol_num(n), k_num(n), g(50,n)

      real*4     PI
      real*4     r, vbin, nid, bin, volume
      real*4     nrho(50), nr(50)
      real*4     gr(50,n), cn(50,n)
C
      PI = acos(-1.0d0)

      i_rdf = 0

C     Number density    

              do i = 1, i_type
              do j = 1, i_type

              if (j.eq.i.or.j.gt.i) then
 
             i_rdf = i_rdf + 1

         nr(i_rdf) = 0.0d0
       nrho(i_rdf) = real( (mol_num(j)) / volume )
      k_num(i_rdf) = mol_num(i)

              endif

           enddo
           enddo
           
C     Finalizing calculations of radial distributions

           do k = 1, i_rdf
           do l = 0, nbin

          r     = bin * ( l + 0.5 )
          vbin  = ( ( l + 1 )**3 - l**3) * bin**3
          nid   = ( 4.0d0 / 3.0d0 ) * PI * vbin * nrho(k)

        gr(k,l) = real( g(k,l) / ( nfr * k_num(k) * nid ) )
          nr(k) =   nr(k) + gr(k,l) * nid
        cn(k,l) = cn(k,l) + nr(k)

           enddo
           enddo

      end

C     ------------------------------------------------------------------ 
C                                                          Print results

      subroutine print_rdf(bin, nbin, i_rdf, gr, cn, k_unit)

      integer n
      parameter (n=1000000)

      integer i, j, nbin, i_rdf, k_unit
      real*4  bin, r, gr(50,n), cn(50,n)

       write(k_unit,*) '#r       ',('  gr',i,'   cn',i, i = 1, i_rdf)

          do j = 0, nbin
             r = bin * ( j + 0.5 )
       write(k_unit,*) r, (gr(i,j), cn(i,j), i = 1, i_rdf)
       enddo 

      end

C     ------------------------------------------------------------------ 
C                                           Assign proper names to atoms

      subroutine assign_name(atm_name, atm_mass)

      integer      i

      real*4       atm_mass
      character*4  atm_name

      real*4       atom_mass(120)
      character*4  atom_name(120)

      character*1  name_1
      character*2  name_2
      character*3  name_3
      character*4  name_4


C     -----------------------------------------------------------------
C          Chemical elements of the periodic table, atomic mass (a.u.). 
C

      atom_mass(1)    =   1.0079         ! H   
      atom_mass(2)    =   4.0026         ! He
      atom_mass(3)    =   6.9410         ! Li
      atom_mass(4)    =   9.0122         ! Be
      atom_mass(5)    =  10.8110         ! B
      atom_mass(6)    =  12.0107         ! C
      atom_mass(7)    =  14.0067         ! N
      atom_mass(8)    =  15.9994         ! O
      atom_mass(9)    =  18.9984         ! F
      atom_mass(10)   =  20.1797         ! Ne
      atom_mass(11)   =  22.9897         ! Na
      atom_mass(12)   =  24.3050         ! Mg
      atom_mass(13)   =  26.9815         ! Al
      atom_mass(14)   =  28.0855         ! Si
      atom_mass(15)   =  30.9738         ! P
      atom_mass(16)   =  32.0650         ! S
      atom_mass(17)   =  35.4530         ! Cl
      atom_mass(18)   =  39.0983         ! K
      atom_mass(19)   =  39.9480         ! Ar
      atom_mass(20)   =  40.0780         ! Ca

      atom_mass(21)   =  44.9559         ! Sc  
      atom_mass(22)   =  47.8670         ! Ti
      atom_mass(23)   =  50.9415         ! V
      atom_mass(24)   =  51.9961         ! Cr 
      atom_mass(25)   =  54.9380         ! Mn 
      atom_mass(26)   =  55.8450         ! Fe 
      atom_mass(27)   =  58.6934         ! Ni 
      atom_mass(28)   =  58.9332         ! Co 
      atom_mass(29)   =  63.5460         ! Cu 
      atom_mass(30)   =  65.3900         ! Zn 
      atom_mass(31)   =  69.7230         ! Ga 
      atom_mass(32)   =  72.6400         ! Ge 
      atom_mass(33)   =  74.9216         ! As 
      atom_mass(34)   =  78.9600         ! Se 
      atom_mass(35)   =  79.9040         ! Br 
      atom_mass(36)   =  83.8000         ! Kr 
      atom_mass(37)   =  85.4678         ! Ru 
      atom_mass(38)   =  87.6200         ! St 
      atom_mass(39)   =  88.9059         ! Y 
      atom_mass(40)   =  91.2240         ! Zr 

      atom_mass(41)   =  92.9064         ! Nb   
      atom_mass(42)   =  95.9400         ! Mo
      atom_mass(43)   =  98.0000         ! Tc
      atom_mass(44)   =  101.070         ! Ru
      atom_mass(45)   =  102.9055        ! Rh
      atom_mass(46)   =  106.4200        ! Pd
      atom_mass(47)   =  107.8682        ! Ag
      atom_mass(48)   =  112.4110        ! Cd
      atom_mass(49)   =  114.8180        ! In
      atom_mass(50)   =  118.7100        ! Sn
      atom_mass(51)   =  121.7600        ! Sb
      atom_mass(52)   =  126.9045        ! I
      atom_mass(53)   =  127.6000        ! Te
      atom_mass(54)   =  131.2930        ! Xe
      atom_mass(55)   =  132.9055        ! Cs
      atom_mass(56)   =  137.3270        ! Ba
      atom_mass(57)   =  138.9055        ! La
      atom_mass(58)   =  140.1160        ! Ce
      atom_mass(59)   =  140.9077        ! Pr
      atom_mass(60)   =  144.2400        ! Nd

      atom_mass(61)   =  145.0000        ! Pm
      atom_mass(62)   =  150.3600        ! Sm
      atom_mass(63)   =  151.9640        ! Eu
      atom_mass(64)   =  157.2500        ! Gd
      atom_mass(65)   =  158.9253        ! Tb
      atom_mass(66)   =  162.5000        ! Dy
      atom_mass(67)   =  164.9303        ! Ho
      atom_mass(68)   =  167.2590        ! Er
      atom_mass(69)   =  168.9342        ! Tm
      atom_mass(70)   =  173.0400        ! Yb
      atom_mass(71)   =  174.9670        ! Lu
      atom_mass(72)   =  178.4900        ! Hf
      atom_mass(73)   =  180.9479        ! Ta
      atom_mass(74)   =  183.8400        ! W
      atom_mass(75)   =  186.2070        ! Re
      atom_mass(76)   =  190.2300        ! Os
      atom_mass(77)   =  192.2170        ! Ir
      atom_mass(78)   =  195.0780        ! Pt
      atom_mass(79)   =  196.9665        ! Au
      atom_mass(80)   =  200.5900        ! Hg

      atom_mass(81)   =  204.3833        ! Ti
      atom_mass(82)   =  207.2000        ! Pb
      atom_mass(83)   =  208.9804        ! Bi
      atom_mass(84)   =  209.0000        ! Po
      atom_mass(85)   =  210.0000        ! At
      atom_mass(86)   =  222.0000        ! Rn
      atom_mass(87)   =  223.0000        ! Fr
      atom_mass(88)   =  226.0000        ! Ra
      atom_mass(89)   =  227.0000        ! Ac
      atom_mass(90)   =  231.0359        ! Pa
      atom_mass(91)   =  232.0381        ! Th
      atom_mass(92)   =  237.0000        ! Np
      atom_mass(93)   =  238.0289        ! U
      atom_mass(94)   =  243.0000        ! Am
      atom_mass(95)   =  244.0000        ! Pu
      atom_mass(96)   =  247.0000        ! Cm
      atom_mass(97)   =  247.0000        ! Bk
      atom_mass(98)   =  251.0000        ! Cf
      atom_mass(99)   =  252.0000        ! Es
      atom_mass(100)  =  257.0000        ! Fm

      atom_mass(101)  =  258.0000        ! Md
      atom_mass(102)  =  259.0000        ! No
      atom_mass(103)  =  261.0000        ! Rf
      atom_mass(104)  =  262.0000        ! Lr
      atom_mass(105)  =  262.0000        ! Db
      atom_mass(106)  =  264.0000        ! Bh
      atom_mass(107)  =  266.0000        ! Sg
      atom_mass(108)  =  268.0000        ! Mt
      atom_mass(109)  =  272.0000        ! Rg
      atom_mass(110)  =  277.0000        ! Hs
      atom_mass(111)  =    0.0000        ! Ds
      atom_mass(112)  =    0.0000        ! Uub
      atom_mass(113)  =    0.0000        ! Uut
      atom_mass(114)  =    0.0000        ! Uuq
      atom_mass(115)  =    0.0000        ! Uup
      atom_mass(116)  =    0.0000        ! Uuh
      atom_mass(117)  =    0.0000        ! Uus
      atom_mass(118)  =    0.0000        ! Uuo
      atom_mass(119)  =    0.0000        ! ---
      atom_mass(120)  =    0.0000        ! ---

C     -----------------------------------------------------------------
C                   Chemical elements of the periodic table, atom names 

      atom_name(1)    =  'H'
      atom_name(2)    =  'He'
      atom_name(3)    =  'Li'
      atom_name(4)    =  'Be'
      atom_name(5)    =  'B'
      atom_name(6)    =  'C'
      atom_name(7)    =  'N'
      atom_name(8)    =  'O'
      atom_name(9)    =  'F'
      atom_name(10)   =  'Ne'
      atom_name(11)   =  'Na'
      atom_name(12)   =  'Mg'
      atom_name(13)   =  'Al'
      atom_name(14)   =  'Si'
      atom_name(15)   =  'P'
      atom_name(16)   =  'S'
      atom_name(17)   =  'Cl'
      atom_name(18)   =  'K'
      atom_name(19)   =  'Ar'
      atom_name(20)   =  'Ca'

      atom_name(21)   =  'Sc'
      atom_name(22)   =  'Ti'
      atom_name(23)   =  'V'
      atom_name(24)   =  'Cr'
      atom_name(25)   =  'Mn'
      atom_name(26)   =  'Fe'
      atom_name(27)   =  'Ni'
      atom_name(28)   =  'Co'
      atom_name(29)   =  'Cu'
      atom_name(30)   =  'Zn'
      atom_name(31)   =  'Ga'
      atom_name(32)   =  'Ge'
      atom_name(33)   =  'As'
      atom_name(34)   =  'Se'
      atom_name(35)   =  'Br'
      atom_name(36)   =  'Kr'
      atom_name(37)   =  'Ru'
      atom_name(38)   =  'St'
      atom_name(39)   =  'Y'
      atom_name(40)   =  'Zr'

      atom_name(41)   =  'Nb'
      atom_name(42)   =  'Mo'
      atom_name(43)   =  'Tc'
      atom_name(44)   =  'Ru'
      atom_name(45)   =  'Rh'
      atom_name(46)   =  'Pd'
      atom_name(47)   =  'Ag'
      atom_name(48)   =  'Cd'
      atom_name(49)   =  'In'
      atom_name(50)   =  'Sn'
      atom_name(51)   =  'Sb'
      atom_name(52)   =  'I'
      atom_name(53)   =  'Te'
      atom_name(54)   =  'Xe'
      atom_name(55)   =  'Cs'
      atom_name(56)   =  'Ba'
      atom_name(57)   =  'La'
      atom_name(58)   =  'Ce'
      atom_name(59)   =  'Pr'
      atom_name(60)   =  'Nd'

      atom_name(61)   =  'Pm'
      atom_name(62)   =  'Sm'
      atom_name(63)   =  'Eu'
      atom_name(64)   =  'Gd'
      atom_name(65)   =  'Tb'
      atom_name(66)   =  'Dy'
      atom_name(67)   =  'Ho'
      atom_name(68)   =  'Er'
      atom_name(69)   =  'Tm'
      atom_name(70)   =  'Yb'
      atom_name(71)   =  'Lu'
      atom_name(72)   =  'Hf'
      atom_name(73)   =  'Ta'
      atom_name(74)   =  'W'
      atom_name(75)   =  'Re'
      atom_name(76)   =  'Os'
      atom_name(77)   =  'Ir'
      atom_name(78)   =  'Pt'
      atom_name(79)   =  'Au'
      atom_name(80)   =  'Hg'

      atom_name(81)   =  'Ti'
      atom_name(82)   =  'Pb'
      atom_name(83)   =  'Bi'
      atom_name(84)   =  'Po'
      atom_name(85)   =  'At'
      atom_name(86)   =  'Rn'
      atom_name(87)   =  'Fr'
      atom_name(88)   =  'Ra'
      atom_name(89)   =  'Ac'
      atom_name(90)   =  'Pa'
      atom_name(91)   =  'Th'
      atom_name(92)   =  'Np'
      atom_name(93)   =  'U'
      atom_name(94)   =  'Am'
      atom_name(95)   =  'Pu'
      atom_name(96)   =  'Cm'
      atom_name(97)   =  'Bk'
      atom_name(98)   =  'Cf'
      atom_name(99)   =  'Es'
      atom_name(100)  =  'Fm'

      atom_name(101)  =  'Md'
      atom_name(102)  =  'No'
      atom_name(103)  =  'Rf'
      atom_name(104)  =  'Lr'
      atom_name(105)  =  'Db'
      atom_name(106)  =  'Bh'
      atom_name(107)  =  'Sg'
      atom_name(108)  =  'Mt'
      atom_name(109)  =  'Rg'
      atom_name(110)  =  'Hs'
      atom_name(111)  =  'Ds'
      atom_name(112)  =  'Uub'
      atom_name(113)  =  'Uut'
      atom_name(114)  =  'Uuq'
      atom_name(115)  =  'Uup'
      atom_name(116)  =  'Uuh'
      atom_name(117)  =  'Uus'
      atom_name(118)  =  'Uuo'
      atom_name(119)  =  'NAY'
      atom_name(120)  =  'NAY'

      name_1 = atm_name
      name_2 = atm_name
      name_3 = atm_name
      name_4 = atm_name

         do i = 1, 120

            if ( name_1 == atom_name(i) )  atm_name = atom_name(i)
            if ( name_2 == atom_name(i) )  atm_name = atom_name(i)
            if ( name_3 == atom_name(i) )  atm_name = atom_name(i)
            if ( name_4 == atom_name(i) )  atm_name = atom_name(i)

            if ( name_1 == atom_name(i) )  atm_mass = atom_mass(i)
            if ( name_2 == atom_name(i) )  atm_mass = atom_mass(i)
            if ( name_3 == atom_name(i) )  atm_mass = atom_mass(i)
            if ( name_4 == atom_name(i) )  atm_mass = atom_mass(i)

         enddo

      end
