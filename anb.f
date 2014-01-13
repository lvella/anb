      program anb
c
************************************************************************
*                                                                      *
*     anb                                                              *
*                                ... means: anb's not best             *
*                                           -     -   -                *
*     A lattice Boltzmann CFD teaching code.                           *
*                                                                      *
*             Joerg BERNSDORF                                          *
*             Contact: research AT bernsdorf.org                       *
*                                                                      *
*             WWW: http://bernsdorf.org/research                       *
*                                                                      *
*     Corrections and hints from Thomas Zeiser and Wolfgang Hamm       *
*     are gratefully acknowledged.                                     *
*                                                                      *
*     Copyright (C) 1998-2008  Joerg BERNSDORF                         *
*                                                                      *
*     This program is free software; you can redistribute it and/or    *
*     modify it under the terms of the GNU General Public License as   *
*     published by the Free Software Foundation; either version 2 of   *
*     the License, or (at your option) any later version.              *
*                                                                      *
*     This program is distributed in the hope that it will be useful,  *
*     but WITHOUT ANY WARRANTY; without even the implied warranty of   *
*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
*     GNU General Public License for more details.                     *
*                                                                      *
*     You should have received a copy of the GNU General Public        *
*     License along with this program; if not, write to the Free       *
*     Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,    *
*     USA.                                                             *
*                                                                      *
*                                                                      *
*     Last change: 2008/05/16                                          *
*                                                                      *
************************************************************************
c
c.....short introduction
c
c.....anb simulates incompressible viscous flow as governed by the 
c     Navier-Stokes equations.
c     But anb is not a Navier-Stokes solver. 
c     anb is a lattice Boltzmann solver, which means, the velocity 
c     discrete boltzmann equation is solved here, using a single time
c     relaxation (BGK) collision operator.
c
c     For more details on this method read:
c             
c                         Y.H.Quian, D.D'Humieres and P.Lallemand,
c                         Lattice BGK Models for Navier-Stokes Equation
c                         Europhys. Lett., 17 (6), pp. 479-484 (1992)
c
c     Other methods belonging to the lattice gas / lattice Boltzmann
c     familiy exist, but this one is very simple, quite good and used
c     by lots of scientists working on this topic. 
c
c     This is a very simple implementation of the lattice BGK scheme,
c     not a very efficient and not a very memory safing one. 
c     *  Do not misunderstand this as a good proposal for writing an 
c        efficient lattice Boltzmann code ! 
c     *  Do not use this code for memory- or time-intensive 
c        computations, it is not even a research code !
c     *  Don't even think about commercial application ! 
c
c     This code was written to show beginners in a simple and 
c     short way the relevant procedures of a lattice Boltzmann solver,
c     pointing on how everything works "in principle". Nearly all 
c     procedures could be implemented other (and better) as it is done
c     here, and even the algorithms used here could be changed to 
c     save memory and increase performance. But the code works correct,
c     and we hope it will be good starting point for the first steps
c     in the lattice Boltzmann field. Good luck ! 
c
      implicit none
c     
c.....parameters
c
c.....grid size in x- and y-dimension
c
      integer  lx,ly
c
c.....ENTER APPROPRIATE VALUES HERE, TAKE CARE: SIZE DOES MATTER !!
c
      parameter(lx=30,ly=20)
c
c.....variables
c
c.....fluid density per link
c
      real*8  density
c
c.....relaxation parameter
c
      real*8  omega 
c
c.....accelleration 
c
      real*8  accel
c
c.....maximum number of iterations
c
      integer  t_max
c
c.....linear dimension for Reynolds number
c
      real*8 r_rey
c
c.....iteration counter
c
      integer time
c
c.....error flag
c
      logical  error
c
c.....obstacle array
c
      logical  obst(lx,ly)
c
c.....fluid densities
c     a 9-speed lattice is used here, other geometries are possible
c
c     the densities are numbered as follows:
c
c              6   2   5
c                \ | /
c              3 - 0 - 1
c                / | \
c              7   4   8
c
c    the lattice nodes are numbered as follows:
c
c
c     ^
c     |
c     y
c
c        :    :    :
c  
c   3    *    *    *  ..
c      
c   2    *    *    *  ..
c      
c   1    *    *    *  ..
c                         x ->
c        1    2    3 
c
c
      real*8  node(0:8,lx,ly)
c
c.....help array for temporarely storage of fluid densities
c
      real*8  n_hlp(0:8,lx,ly)
c
c.....average velocity, computed by subroutine 'write_velocity'
c
      real*8  vel
c
c.....startup information message     
c
      write (6,*)
      write (6,*) 'anb (2008/05/16)'
      write (6,*) 'Copyright (C) 1998-2008 Joerg Bernsdorf /'
      write (6,*) 'anb comes with ABSOLUTELY NO WARRANTY;' 
      write (6,*) 'This is free software, and you are welcome to redistr
     &ibute it' 
      write (6,*) 'under certain conditions; see the file COPYING for de
     &tails.'
      write (6,*) 'All rights reserved.'
c
      write (6,*)
      write (6,*) '****************************************************'
      write (6,*) '***              anb starting ...                ***'
      write (6,*) '****************************************************'
      write (6,*) '*** Precompiled for lattice size lx = ',lx 
      write (6,*) '***                              ly = ',ly 
      write (6,*) '****************************************************'
      write (6,*) '***'
c
c=======================================================================
c     begin initialisation
c=======================================================================
c
c.....initialize error flag
c
      error = .false.
c
c.....read parameter file
c.....in this file you can enter all relevant parameters,  
c     only lattice size must be fixed befor compilation !
c
      call read_parametrs(error,t_max,density,accel,omega,r_rey)
c
c.....if an I/O error occurs while reading the parameter file, 
c     the "error"-flag is set "true" and the program stops.
c
      if (error) goto 990
c
c.....read obstacle file
c.....in this file you can enter the x-,and y-coordinates of 
c     any obstacles,wall boundaries are also defined here by 
c     adding single obstacles.
c
      call read_obstacles(error,obst,lx,ly)     
c
c.....if an I/O error occurs reading while the obstacle file, 
c     the "error"-flag is set "true" and the program stops.
c
      if (error) goto 990
c
      call init_density(lx,ly,density,node)
c
c.....mean sqaure flow is monitored and stored in the file'anb_qxm.out',
c     one value for each iteration. The file is opened here.
c
      open(10,file='anb_qxm.out')      
c
c=======================================================================
c     end initialisation
c=======================================================================
c=======================================================================
c     begin iterations
c=======================================================================
c
c.....main loop
c
      do 100 time = 1, t_max
c
c.......the integral fluid density is checke each t_max/10 iteration.
c       this is a good indicator, if the program is going to crash ...
c       The integral fluid density should be constant all time ...
c
        if (time .ge. 10 .and. mod(time,t_max/10) .eq. 0) then
c        
          call check_density(lx,ly,node,time)
c
        end if
c
c.......directed flow is induced by density redistribution in the first 
c       lattice column. This is not too clever, since the resulting
c       reynolds number can not be controlled and reaching stady state 
c       takes quite some time, but it is simple and it works ...
c
        call redistribute(lx,ly,obst,node,accel,density)      
c
c.......density propagation: all fluid densities are propagated from
c       non-occupied nodes along the lattice connection lines
c       to their next neighbours, periodic boundary conditiones are 
c       applied in each direction.
c
        call propagate(lx,ly,node,n_hlp)
c
c.......bounc back from obstacles: this is the no-slip boundary-
c       condition.
c       The velocity vector of all fluid densities is inverted, so all
c       the fluid densities will be sent back to the node  where they 
c       were located before the last propagation step, but with opposite
c       velocity vector 
c       ... there exist lots of other possibilities.
c
        call bounceback(lx,ly,obst,node,n_hlp)
c
c.......density relaxation: a single time relaxation with relaxation 
c       parameter omega is applied here. This step is only "local", 
c       nothing is propagated through the lattice.
c
        call relaxation(density,omega,lx,ly,node,n_hlp,obst)
c
c.......average flow velocity is computed at a cross section in the  
c       middle of the channel and written to the file "anb_qx.out"
c       ever iteration. this is a good controll for the convergence of
c       the program.
c
        call write_velocity(lx,ly,time,obst,node,vel)        
c
c.....end of the main loop
c
  100 continue
c
c.....compute fluid-velocities u,v and pressure from velocity 
c     distribution, and write to file anb_rs.out.
c
      call write_results(lx,ly,obst,node,density)
c
c.....compute reynolds number
c
      call comp_rey(lx,ly,obst,node,time,omega,density,r_rey)
c
      goto 999
c
c.....here we get only, if the "error" flag was set "true",  something
c     went wrong reading the files .The program stops with an error
c     message.
c
  990 write (6,*) '!!! error: program stopped during iteration =', time
      write (6,*) '!!!'
c      
  999 continue
c
c.....the file for mean square flow output (anb_qx.out) is closed. 
c     Look at this file before starting any other evaluating, it tells 
c     you, if you reached stady state and if you got reasonable results!
c
      close(10)
c
      write (6,*) '********************    end     ********************'
c
c      
      stop
      end	
c
c
      subroutine read_parametrs(error,t_max,density,accel,omega,r_rey)
c
************************************************************************
*                                                                      *
*     Input run-time parameters from file 'anb.par'                    *
*                                                                      *
*             Joerg BERNSDORF                                          *
*             C&C Research Laboratories, NEC Europe Ltd.               *
*             Rathausalle 10                                           *
*             D-53757 Sankt Augustin, Germany                          *
*                                                                      *
*     Last change: 1999/09/28                                          *
*                                                                      *
************************************************************************
c
      implicit none
c
      real*8  density,accel,omega,r_rey
c
      integer  t_max
c
      logical  error
c
c.....open parameter file
c
      open(10,file='anb.par')
c
c.......line 1: number of iterations
c
        read(10,*,err=900) t_max
c
c.......line 2: fluid density per link
c
        read(10,*,err=900) density
c
c.......line 3: density redistribution
c
        read(10,*,err=900) accel
c
c.......line 4: relaxation parameter
c
        read(10,*,err=900) omega
c
c.......line 5: linear dimension (for reynolds number)
c
        read(10,*,err=900) r_rey
c
c.....close parameter file
c
      close(10)      
c
c.....information message
c
      write (6,*) '*** Paramters read from file anb.par.'
      write (6,*) '***'
c
      goto 999
c
c.....error message: file read error
c
  900 write (6,*) '!!! Error reading file anb.par'
      write (6,*) '!!!'
c
      goto 990
c
  990 error = .true.
c
  999 continue
c
c
      return
      end
c
c
      subroutine read_obstacles(error,obst,lx,ly)
c
************************************************************************
*                                                                      *
*     Input obstacle file 'anb.obs'                                    *
*                                                                      *
*             Joerg BERNSDORF                                          *
*             C&C Research Laboratories, NEC Europe Ltd.               *
*             Rathausalle 10                                           *
*             D-53757 Sankt Augustin, Germany                          *
*                                                                      *
*     Last change: 1998/08/25                                          *
*                                                                      *
************************************************************************
c
      implicit none
c     
      integer  lx,ly
c
      logical  error,obst(lx,ly)
c
c.....local variables
c
      integer  x,y
c
c.....no obstacles in obstacle array
c
        do 10 y = 1, ly
          do 10 x = 1, lx
c
   10       obst(x,y) = .false.
c
c.....open obstacle file
c.....the obstacle file is defined by the x- and y- coordinates of the
c     obstacles. Each obstacle has a line of its own. Also boundaries
c     have to be defined here.
c
      open(10,file='anb.obs')
c
c.....read obstacle coordinates
c
   20 continue
c
        read(10,*,end=50,err=900) x,y
c
c.......check if obstacle inside domain boundaries
c
        if (x .le. lx .and. y .le. ly) then
c
c.......define obstacle
c
          obst(x,y) = .true.
c
        else
c
          write(6,*) '!!! Obstacle out of range, skipped'
          write(6,*) '!!! lx = ', x, ' , ly = ', y
          write(6,*) '!!!'
c
        end if
c
      goto 20
c
   50 continue
c
c.....close obstacle file
c
      close(10)
c
      write (6,*) '*** Geometry information read from file anb.obs.'
      write (6,*) '***'
c
      goto 999
c
c.....error message: file read error
c
  900 write (6,*) '!!! Error reading file anb.obs'
      write (6,*) ' '
c
      goto 990
c
  990 error = .true.
c
  999 continue
c
c
      return
      end
c
c
      subroutine init_density(lx,ly,density,node)
c
************************************************************************
*                                                                      *
*     Initialize density distribution function n with equilibrium      *
*     for zero velocity                                                *
*                                                                      *
*             Joerg BERNSDORF                                          *
*             C&C Research Laboratories, NEC Europe Ltd.               *
*             Rathausalle 10                                           *
*             D-53757 Sankt Augustin, Germany                          *
*                                                                      *
*     Last change: 1999/09/28                                          *
*                                                                      *
************************************************************************
c
      implicit none
c
      integer  lx,ly
c
      real*8  density,node(0:8,lx,ly)
c
c.....local variables
c
      integer  x,y
      real*8  t_0,t_1,t_2
c
c.....compute weighting factors (depending on lattice geometry)
c
      t_0 = density *  4.d0 / 9.d0
      t_1 = density /  9.d0
      t_2 = density / 36.d0
c
c.....loop over computational domain
c
      do 10 x = 1, lx
        do 10 y = 1, ly
c
c.........zero velocity density
c
          node(0,x,y) = t_0
c
c.........equilibrium densities for axis speeds
c
          node(1,x,y) = t_1
          node(2,x,y) = t_1
          node(3,x,y) = t_1
          node(4,x,y) = t_1
c
c.........equilibrium densities for diagonal speeds
c
          node(5,x,y) = t_2
          node(6,x,y) = t_2
          node(7,x,y) = t_2
          node(8,x,y) = t_2
c
   10 continue
c
c
      return
      end
c
c
      subroutine check_density(lx,ly,node,time)  
c
************************************************************************
*                                                                      *
*     compute integral density                                         *
*                                                                      *
*             Joerg BERNSDORF                                          *
*             C&C Research Laboratories, NEC Europe Ltd.               *
*             Rathausalle 10                                           *
*             D-53757 Sankt Augustin, Germany                          *
*                                                                      *
*     Last change: 2000/01/14                                          *
*                                                                      *
************************************************************************
c
      implicit none
c     
      integer  lx,ly,time
c
      real*8  node(0:8,lx,ly)
c
c
c.....local variables
c
      integer  x,y,n
c
      real*8 n_sum
c
c
      n_sum = 0.d0
c
c.....loop over computational domain
c
        do 10 y = 1, ly
          do 10 x = 1, lx
c
c...........loop over all densities
c
            do 10 n = 0, 8
c
c.............sum up densities
c
   10         n_sum = n_sum + node(n,x,y)
c
      write(6,*) '*** Iteration number = ', time
      write(6,*) '*** Integral density = ', n_sum
      write(6,*) '***'
c
c
      return
      end
c
c
      subroutine redistribute(lx,ly,obst,node,accel,density)
c
************************************************************************
*                                                                      *
*     density redistribution in first lattice column                   *
*                                                                      *
*             Joerg BERNSDORF                                          *
*             C&C Research Laboratories, NEC Europe Ltd.               *
*             Rathausalle 10                                           *
*             D-53757 Sankt Augustin, Germany                          *
*                                                                      *
*     Last change: 1999/09/28                                          *
*                                                                      *
************************************************************************
c
      implicit none
c
      integer  lx,ly
c
      logical  obst(lx,ly)
c
      real*8   node(0:8,lx,ly),accel,density
c
c.....local variables
c
      integer  y
      real*8  t_1,t_2
c
c.....compute weighting factors (depending on lattice geometry) for 
c     increasing/decreasing inlet densities
c
      t_1 = density * accel /  9.d0
      t_2 = density * accel / 36.d0

      do 10 y = 1, ly
c
c.......accelerate flow only on non-occupied nodes
c
        if (.not. obst(1,y) .and.
c
c.........check to avoid negative densities
c
     &    node(3,1,y) - t_1 .gt. 0. .and.
     &    node(6,1,y) - t_2 .gt. 0. .and.
     &    node(7,1,y) - t_2 .gt. 0.) then 
c
c.........increase east
c
          node(1,1,y) = node(1,1,y) + t_1
c
c.........decrease west
c
           node(3,1,y) = node(3,1,y) - t_1
c
c.........increase north-east
c
          node(5,1,y) = node(5,1,y) + t_2
c
c.........decrease north-west
c
          node(6,1,y) = node(6,1,y) - t_2
c
c.........decrease south-west
c
          node(7,1,y) = node(7,1,y) - t_2
c
c.........increase south-east
c
          node(8,1,y) = node(8,1,y) + t_2
c
        end if
c
   10 continue
c
c
      return
      end
c
c
      subroutine bounceback(lx,ly,obst,node,n_hlp)
c
************************************************************************
*                                                                      *
*     Fluid densities are rotated. By the next propagation step, this  *
*     results in a bounce back from obstacle nodes.                    *
*                                                                      *
*             Joerg BERNSDORF                                          *
*             C&C Research Laboratories, NEC Europe Ltd.               *
*             Rathausalle 10                                           *
*             D-53757 Sankt Augustin, Germany                          *
*                                                                      *
*     Last change: 1999/09/28                                          *
*                                                                      *
************************************************************************
c
      implicit none
c
      integer  lx,ly
c
      logical  obst(lx,ly)
c
      real*8   node(0:8,lx,ly),n_hlp(0:8,lx,ly)
c
c.....local variables
c
      integer  x,y
c
c.....loop over all nodes
c
      do 10 x = 1, lx
        do 10 y = 1, ly
c
c.........consider only obstacle nodes
c
          if (obst(x,y)) then
c
c...........rotate all ensities and write back to node
c
c...........east
c
            node(1,x,y) = n_hlp(3,x,y)
c
c...........north
c
            node(2,x,y) = n_hlp(4,x,y)
c
c...........west
c
            node(3,x,y) = n_hlp(1,x,y)
c
c...........south
c
            node(4,x,y) = n_hlp(2,x,y)
c
c...........north-east
c
            node(5,x,y) = n_hlp(7,x,y)
c
c...........north-west
c
            node(6,x,y) = n_hlp(8,x,y)
c
c...........south-west
c
            node(7,x,y) = n_hlp(5,x,y)
c
c...........south-east
c
            node(8,x,y) = n_hlp(6,x,y)
c
          end if
c
   10 continue
c
c
      return
      end
c
c
      subroutine propagate(lx,ly,node,n_hlp)
c
************************************************************************
*                                                                      *
*     Propagate fluid densities to their next neighbour nodes          *
*                                                                      *
*             Joerg BERNSDORF                                          *
*             C&C Research Laboratories, NEC Europe Ltd.               *
*             Rathausalle 10                                           *
*             D-53757 Sankt Augustin, Germany                          *
*                                                                      *
*     Last change: 1999/09/28                                          *
*                                                                      *
************************************************************************
c
      implicit none
c
      integer  lx,ly
c
      real*8   node(0:8,lx,ly),n_hlp(0:8,lx,ly)
c
c.....local variables
c
      integer  x,y,x_e,x_w,y_n,y_s
c
c.....loop over all nodes
c
      do 10 x = 1, lx
        do 10 y = 1, ly
c
c.........compute upper and right next neighbour nodes with regard
c         to periodic boundaries
c
          y_n = mod(y,ly) + 1
          x_e = mod(x,lx) + 1
c
c.........compute lower and left next neighbour nodes with regard to
c         periodic boundaries

          y_s = ly - mod(ly + 1 - y, ly)
          x_w = lx - mod(lx + 1 - x, lx)
c
c.........density propagation
c
c.........zero: just copy
c
          n_hlp(0,x  ,y  ) = node(0,x,y)
c
c.........east
c
          n_hlp(1,x_e,y  ) = node(1,x,y)
c
c.........north
c
          n_hlp(2,x  ,y_n) = node(2,x,y)
c
c.........west
c
          n_hlp(3,x_w,y  ) = node(3,x,y)
c
c.........south
c
          n_hlp(4,x  ,y_s) = node(4,x,y)
c
c.........north-east
c
          n_hlp(5,x_e,y_n) = node(5,x,y)
c
c.........north-west
c
          n_hlp(6,x_w,y_n) = node(6,x,y)
c
c.........south-west
c
          n_hlp(7,x_w,y_s) = node(7,x,y)
c
c.........south-east
c
          n_hlp(8,x_e,y_s) = node(8,x,y)
c
   10 continue
c
c
      return
      end
c
c
      subroutine relaxation(density,omega,lx,ly,node,n_hlp,obst)
c
************************************************************************
*                                                                      *
*     One-step density relaxation process                              *
*                                                                      *
*             Joerg BERNSDORF                                          *
*             C&C Research Laboratories, NEC Europe Ltd.               *
*             Rathausalle 10                                           *
*             D-53757 Sankt Augustin, Germany                          *
*                                                                      *
*     Last change: 1999/09/28                                          *
*                                                                      *
************************************************************************
c
      implicit none
c
      integer  lx,ly
c
      logical  obst(lx,ly)
c
      real*8   density,omega,node(0:8,lx,ly),n_hlp(0:8,lx,ly)
c
c.....local variables
c
      integer  x,y,i
c
      real*8  c_squ,t_0,t_1,t_2,u_x,u_y,u_n(8),n_equ(0:8),u_squ,d_loc
c
c.....weighting factors (depending on lattice geometry)
c
      t_0 = 4.d0 / 9.d0
      t_1 = 1.d0 /  9.d0
      t_2 = 1.d0 / 36.d0
c
c.....square speed of sound
c
      c_squ = 1.d0 / 3.d0
c
c.....loop over all nodes
c.....attention: actual densities are stored after the propagation
c                step in the help-array n_hlp !
c
      do 10 x = 1, lx
        do 10 y = 1, ly
c
c.........only free nodes are considered here
c
          if (.not. obst(x,y)) then
c
c...........integral local density
c
c...........initialize variable d_loc
c
            d_loc = 0.d0
c
            do 20 i = 0, 8 
c
              d_loc = d_loc + n_hlp(i,x,y)
c
   20       continue
c
c...........x-, and y- velocity components
c
            u_x = (n_hlp(1,x,y) + n_hlp(5,x,y) + n_hlp(8,x,y)
     &           -(n_hlp(3,x,y) + n_hlp(6,x,y) + n_hlp(7,x,y))) / d_loc
c
            u_y = (n_hlp(2,x,y) + n_hlp(5,x,y) + n_hlp(6,x,y)
     &           -(n_hlp(4,x,y) + n_hlp(7,x,y) + n_hlp(8,x,y))) / d_loc
c
c...........square velocity
c
            u_squ = u_x * u_x + u_y * u_y
c
c...........n- velocity compnents (n = lattice node connection vectors)
c...........this is only necessary for clearence, and only 3 speeds would
c...........be necessary
c
            u_n(1) =   u_x
            u_n(2) =         u_y
            u_n(3) = - u_x
            u_n(4) =       - u_y
            u_n(5) =   u_x + u_y
            u_n(6) = - u_x + u_y
            u_n(7) = - u_x - u_y
            u_n(8) =   u_x - u_y
c
c...........equilibrium densities
c...........this can be rewritten to improve computational performance
c...........considerabely !
c
c...........zero velocity density
c
            n_equ(0) = t_0 * d_loc * (1.d0 - u_squ / (2.d0 * c_squ))
c
c...........axis speeds (factor: t_1)
c
            n_equ(1) = t_1 * d_loc * (1.d0 + u_n(1) / c_squ 
     &               + u_n(1) ** 2.d0 / (2.d0 * c_squ ** 2.d0) 
     &               - u_squ / (2.d0 * c_squ))
c
            n_equ(2) = t_1 * d_loc * (1.d0 + u_n(2) / c_squ 
     &               + u_n(2) ** 2.d0 / (2.d0 * c_squ ** 2.d0) 
     &               - u_squ / (2.d0 * c_squ))
c
            n_equ(3) = t_1 * d_loc * (1.d0 + u_n(3) / c_squ 
     &               + u_n(3) ** 2.d0 / (2.d0 * c_squ ** 2.d0) 
     &               - u_squ / (2.d0 * c_squ))
c
            n_equ(4) = t_1 * d_loc * (1.d0 + u_n(4) / c_squ 
     &               + u_n(4) ** 2.d0 / (2.d0 * c_squ ** 2.d0) 
     &               - u_squ / (2.d0 * c_squ))
c
c...........diagonal speeds (factor: t_2)
c
            n_equ(5) = t_2 * d_loc * (1.d0 + u_n(5) / c_squ 
     &               + u_n(5) ** 2.d0 / (2.d0 * c_squ ** 2.d0) 
     &               - u_squ / (2.d0 * c_squ))
c
            n_equ(6) = t_2 * d_loc * (1.d0 + u_n(6) / c_squ 
     &               + u_n(6) ** 2.d0 / (2.d0 * c_squ ** 2.d0) 
     &               - u_squ / (2.d0 * c_squ))
c
            n_equ(7) = t_2 * d_loc * (1.d0 + u_n(7) / c_squ 
     &               + u_n(7) ** 2.d0 / (2.d0 * c_squ ** 2.d0) 
     &               - u_squ / (2.d0 * c_squ))
c
            n_equ(8) = t_2 * d_loc * (1.d0 + u_n(8) / c_squ 
     &               + u_n(8) ** 2.d0 / (2.d0 * c_squ ** 2.d0) 
     &               - u_squ / (2.d0 * c_squ))
c
c...........relaxation step
c
            do 30 i = 0, 8
c
              node(i,x,y) = n_hlp(i,x,y) 
     &                      + omega * (n_equ(i) - n_hlp(i,x,y))
c
   30       continue
c
          end if
c
   10 continue
c
c
      return
      end
c
c
      subroutine write_velocity(lx,ly,time,obst,node,vel)        
c
************************************************************************
*                                                                      *
*     compute average velocity                                         *
*                                                                      *
*             Joerg BERNSDORF                                          *
*             C&C Research Laboratories, NEC Europe Ltd.               *
*             Rathausalle 10                                           *
*             D-53757 Sankt Augustin, Germany                          *
*                                                                      *
*     Last change: 1999/09/28                                          *
*                                                                      *
************************************************************************
c
      implicit none
c
      integer  lx,ly,time
c
      logical  obst(lx,ly)
c
      real*8   node(0:8,lx,ly),vel
c
c.....local variables
c
      integer  x,y,i,n_free
c
      real*8  u_x,d_loc
c
c.....loop over channel cross section at half channel length lx/2
c
      x = int(float(lx) / 2.d0)
c
c.....initialize counter
c
      n_free = 0
      u_x = 0.d0
c
      do 10 y = 1, ly
c
c.......only non-occupied nodes are considered here
c
        if(.not. obst(x,y)) then
c
c.........integral local density
c
c.........initialize variable d_loc
c
          d_loc = 0.d0
c
          do 20 i = 0, 8 
c
            d_loc = d_loc + node(i,x,y)
c
   20     continue
c
c.........x-, and y- velocity components
c
          u_x = u_x + (node(1,x,y) + node(5,x,y) + node(8,x,y)
     &               -(node(3,x,y) + node(6,x,y) + node(7,x,y))) / d_loc
c
c.........increase counter
c
          n_free = n_free + 1
c
        end if
c
   10 continue
c
c.....average velocity
c
      vel = u_x / float(n_free)
c
c.....write to file
c
      write(10,*) time, vel
c
c
      return
      end
c
c
      subroutine  write_results(lx,ly,obst,node,density)
c
************************************************************************
*                                                                      *
*     Output of rsults to file 'anb.dat'                               *
*                                                                      *
*             Joerg BERNSDORF                                          *
*             C&C Research Laboratories, NEC Europe Ltd.               *
*             Rathausalle 10                                           *
*             D-53757 Sankt Augustin, Germany                          *
*                                                                      *
*     Last change: 1999/09/28                                          *
*                                                                      *
************************************************************************
c
      implicit none
c     
      integer  lx,ly
c
      real*8  node(0:8,lx,ly),density
c
      logical  obst(lx,ly)
c
c.....local variables
c
      integer  x,y,i,obsval
c
      real*8  u_x,u_y,d_loc,press,c_squ
c
c.....square speed of sound
c
      c_squ = 1.d0 / 3.d0
c
c.....open results output file
c
      open(11,file='anb.dat')
c
c.....write header for postprocessing with TECPLOT software
c.....uncomment following line, if this header should be printed
c
c
      write(11,*) 'VARIABLES = X, Y, VX, VY, PRESS, OBST' 
      write(11,*) 'ZONE I=', lx, ', J=', ly, ', F=POINT'
c
c.....loop over all nodes
c.....attention: actual densities are stored after the propagation
c                step in the help-array n_hlp !
c
      do 10 y = 1, ly
        do 10 x = 1, lx
c
c.........if obstacle node, nothing is to do ...
c
          if (obst(x,y)) then 
c
c...........obstacle indicator
c
            obsval = 1
c
c...........velocity components = 0
c
            u_x = 0.d0
            u_y = 0.d0
c
c...........pressure = average pressure
c
            press = density * c_squ
c
          else
c
c...........integral local density
c
c...........initialize variable d_loc
c
            d_loc = 0.d0
c
            do 20 i = 0, 8 
c
              d_loc = d_loc + node(i,x,y)
c
   20       continue
c
c...........x-, and y- velocity components
c
            u_x = (node(1,x,y) + node(5,x,y) + node(8,x,y)
     &           -(node(3,x,y) + node(6,x,y) + node(7,x,y))) / d_loc
c
            u_y = (node(2,x,y) + node(5,x,y) + node(6,x,y)
     &           -(node(4,x,y) + node(7,x,y) + node(8,x,y))) / d_loc
c
c...........pressure
c
            press = d_loc * c_squ
c 
            obsval = 0
c 
          end if
c
c.........write results to file
c
          write(11,*) x, y, u_x, u_y, press, obsval
c
   10 continue
c
c.....close file 'anb.dat'
c
      close(11)
c
c
      return
      end
c
c
      subroutine comp_rey(lx,ly,obst,node,time,omega,density,r_rey)
c
************************************************************************
*                                                                      *
*     compute Reynolds number                                          *
*                                                                      *
*             Joerg BERNSDORF                                          *
*             C&C Research Laboratories, NEC Europe Ltd.               *
*             Rathausalle 10                                           *
*             D-53757 Sankt Augustin, Germany                          *
*                                                                      *
*     Last change: 1999/09/28                                          *
*                                                                      *
************************************************************************
c
      implicit none
c
      integer  lx,ly,time
c
      real*8  density,node(0:8,lx,ly),omega,r_rey
c
      logical  obst(lx,ly)
c
c.....local variables
c
      real*8  vel,visc,rey
c
c.....compute average velocity
c
      call write_velocity(lx,ly,time,obst,node,vel)        
c
c.....compute viscosity
c
      visc = 1.d0 / 6.d0 * (2.d0 / omega - 1.d0)
c
c.....compute Reynolds number
c
      rey = vel * r_rey / visc
c
c.....messages
c
       write (6,*) '*** Calculations finished, results:'
       write (6,*) '***'
       write (6,*) '*** viscosity = ', visc
       write (6,*) '*** average velocity = ', vel
       write (6,*) '*** Reynolds number = ', rey
       write (6,*) '***'
       write (6,*) '*** In the file anb.dat, you can find local'
       write (6,*) '*** information about the simulated flow.'
       write (6,*) '***'
       write (6,*) '*** In the file anb_qx.out, you can find the average
     &'
       write (6,*) '*** flow velocity plotted as a function of time.'
       write (6,*) '***'
c
c
      return
      end
