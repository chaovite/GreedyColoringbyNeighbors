! ================================================================ 
!                ~~* Greedy Graph Coloring *~~              
! ================================================================ 
!
! A fortran program that implements the Greedy Graph algorithm
! to color a 2D FEM mesh. This is a fortran version of the matlab code
!
! The program takes in a ELEM.dat file and number of threads 
!  from the command line and write a COLOR.dat file.
!
! To compile: 
!     gfortran GreedyColoringbyNeighbors.F90 -o color_exe
!
! To run:
!     ./color_exe ELEM.dat nthreads COLOR.dat [NEB.dat]
!
! ================================================================ 
!           Written by Chao Liang, Nov 7th, 2021
! ================================================================ 
!

PROGRAM GreedyColoringbyNeighbors
  INTEGER :: i, nthreads
  CHARACTER(len=80)   :: ARG, ELEM_FILE, COLOR_FILE, NEB_FILE
  INTEGER, ALLOCATABLE :: ELEM(:,:), NEB(:,:), Colors(:)
  INTEGER :: NEL, NNODE 
  INTEGER, PARAMETER :: max_neighbor = 50 

  ! Get input parameters from command line arguments
  call ReadCommandLineInput(ELEM_FILE, nthreads, COLOR_FILE, NEB_FILE)
  
  ! Read in the first line of ELEM_FILE to obtain NEL, NNODE
  call Get_ELEM_SIZE(NEL, NNODE, ELEM_FILE)
  allocate(ELEM(NNODE, NEL), NEB(max_neighbor, NEL), Colors(NEL))
  
  ! Read ELEM_FILE to obtain ELEM, the connectivity array
  call READ_ELEM(ELEM, ELEM_FILE)

  ! Find neighbors for each element
  call FindNeighbors(ELEM, NEB) 

  ! Color the FEM mesh using neighbors
  call ColorByNeighbors(NEB, nthreads, Colors, NColor)

  !-------------------- writing output files --------------------- 
  ! write Colors into a file
  open(unit=2, FILE=COLOR_FILE, status='replace', action='write')
  write(2,   '(I6)') Colors  
  close(unit=2)

  
  ! write Neighbors into a file
  if (len(TRIM(NEB_FILE))>0) then
      open(unit=1, FILE=NEB_FILE, status='replace', action='write')
      do i =1, NEL
          do j = 1, max_neighbor
              write(1, '(I6, 1X)', advance='no') NEB(j, i) 
          end do
          write(1, *) ''
      end do
      close(unit=1)
  end if
  !--------------------------------------------------------------- 
  deallocate(ELEM, NEB, Colors)

  contains
      subroutine ReadCommandLineInput(ELEM_FILE, nthreads, COLOR_FILE, NEB_FILE)
          CHARACTER(len=80) :: arg
          CHARACTER(len=80), intent(out) :: ELEM_FILE, COLOR_FILE, NEB_FILE
          INTEGER, intent(out) :: nthreads
          INTEGER :: i

          ! set default parameters
          nthreads   = 2
          ELEM_FILE  = 'ELEM.dat'
          COLOR_FILE = 'COLOR.dat'
          NEB_FILE   = '' ! do not output NEB_FILE

          ! read input from the command lines
          DO i = 1, 4
            CALL get_command_argument(i, arg)
            IF (LEN_TRIM(arg) > 0) then
                SELECT CASE(i)
                CASE(1)
                    ! read in the ELEL.dat filename
                    read(arg,*) ELEM_FILE
                    ELEM_FILE = TRIM(ELEM_FILE)
                CASE(2)
                    read(arg,*) nthreads
                CASE(3)
                    ! read in the COLOR.dat filename
                    read(arg,*) COLOR_FILE
                    COLOR_FILE = TRIM(COLOR_FILE)
                CASE(4)
                    ! read in the COLOR.dat filename
                    read(arg,*) NEB_FILE
                    NEB_FILE = TRIM(NEB_FILE)
                END SELECT
            END IF
          END DO
      write(*, '(/)')
      write(*, '(A)')     "======================================="
      write(*, '(A)')     "*                                     *"
      write(*, '(A)')     "*        Greedy Graph Coloring        *"
      write(*, '(A)')     "*                                     *"
      write(*, '(A)')     "======================================="
      write(*, *)
      write(*, '(A)')     "---------------------------------------"
      write(*, '(A)')     "Reading inputs from command line ......"
      write(*, '(A, A)')  "ELEM_FILE  = ", ELEM_FILE
      write(*, '(A, I3)') "NTHREADS   = ", NTHREADS
      write(*, '(A, A)')  "COLOR_FILE = ", COLOR_FILE
      write(*, '(A, A)')  "NEB_FILE   = ", NEB_FILE
      write(*, '(A)')     "---------------------------------------"
      write(*, *)
      end subroutine ReadCommandLineInput

      subroutine READ_ELEM(ELEM, ELEM_FILE)
          INTEGER, DIMENSION(:,:), INTENT(INOUT) :: ELEM
          INTEGER :: NEL, NNODE, i
          CHARACTER(80) :: ELEM_FILE

          NEL   = size(ELEM, 2)
          NNODE = size(ELEM, 1)
          open(UNIT=1, FILE=ELEM_FILE, STATUS='OLD',ACTION='READ')
          read(1, *) NEL, NNODE
          read(1, *) ELEM 
          close(1)
          write(*, '(A)')     "Finish reading ELEM, connectivity array"
          write(*, '(A)')     "---------------------------------------"
      end subroutine READ_ELEM

      ! Get the size of ELEM array from the first line of ELEM_FILE
      subroutine Get_ELEM_SIZE(NEL, NNODE, ELEM_FILE)
          INTEGER, INTENT(OUT):: NEL, NNODE
          CHARACTER(80) :: ELEM_FILE
          open(UNIT=1, FILE=ELEM_FILE, STATUS='OLD',ACTION='READ')
          read(1, *) NEL, NNODE
          write(*, '(A, I8, A, I1)') "NEL = ", NEL, ' NNODE = ', NNODE
          close(UNIT=1)
      end subroutine Get_ELEM_SIZE

      ! find the neighboring elements for each element
      subroutine FindNeighbors(ELEM, NEB)  
          INTEGER, INTENT(IN):: ELEM(:, :)
          INTEGER, INTENT(OUT):: NEB(max_neighbor, size(ELEM, 2))
          ! the first dimension of NEB must be large than 
          ! maximum number of neighbors, I set it to be max_neighbor
          INTEGER :: i, j, k, NEL, i_node, cnt_i, n2e(max_neighbor)
          INTEGER :: e_i, e_j
          INTEGER, allocatable :: node2el(:,:), cnt1(:)
          Integer, allocatable :: tmp(:)
          Integer :: cnt2(size(ELEM, 2)), nodes(size(ELEM, 1))

          NEL     = size(ELEM, 2)
          nnode   = maxval(ELEM) ! maximum node index  

          allocate(node2el(max_neighbor, nnode), cnt1(nnode)) 
          node2el = 0
          cnt1    = 0
          NEB     = -1
          cnt2    = 0

          ! build the array node2el
          do i = 1, NEL
              nodes = ELEM(:, i)
              do j = 1, size(nodes)
                  i_node = nodes(j)
                  ! (i, i_node) should be included in node2el
                  ! if not included then update both count and nodel2el
                  if (allocated(tmp)) deallocate(tmp)
                  tmp  = findloc(node2el(:, i_node),i)
                  if (size(tmp)==1 .and. tmp(1)==0) then
                      ! update the count of elements for each node 
                      cnt1(i_node) = cnt1(i_node) + 1
                      node2el(cnt1(i_node),i_node) = i
                  end if
              end do
          end do

          ! find the neighboring elements
          do i = 1, nnode
              n2e = node2el(:,i)
              cnt_i = cnt1(i) !# elements for each node
              do k = 1, cnt_i
                  e_i = n2e(k)
                  do j = k+1, cnt_i
                      e_j = n2e(j)
                      ! check if e_i, e_j has already been included
                      if (allocated(tmp)) deallocate(tmp)
                      tmp  = findloc(NEB(:, e_i), e_j)
                      if (size(tmp)==1 .and. tmp(1)==0) then
                          ! if not, include pair
                          cnt2(e_i) = cnt2(e_i) + 1
                          cnt2(e_j) = cnt2(e_j) + 1
                          NEB(cnt2(e_i), e_i) = e_j
                          NEB(cnt2(e_j), e_j) = e_i
                      end if
                  end do
              end do
          end do
          deallocate(node2el, cnt1, tmp) 
          write(*,*) "done building neighbors..."
      end subroutine FINDNeighbors

      ! color the mesh for nthreads, return color and number of colors
      subroutine ColorByNeighbors(NEB, nthreads, Colors, NColor)
          INTEGER, INTENT(IN)  :: NEB(:,:), nthreads
          INTEGER, INTENT(OUT) :: NColor, Colors(size(NEB, 2))
          INTEGER :: NEL, i, j, ColorCount(size(NEB, 2))
          INTEGER :: BlockedColors(max_neighbor), Neighbors(max_neighbor)
          LOGICAL :: isfree
          INTEGER, allocatable :: tmp(:)

          NEL        = size(NEB, 2)
          Colors     = 0
          Colors(1)  = 1
          NColor     = 1
          ColorCount = 0
          ColorCount(1) = 1
          
          ! loop through the second to the last element
          do i = 2, NEL
              BlockedColors = -1
              Neighbors = NEB(:, i)
              do j = 1, max_neighbor
                  if (Neighbors(j)<=0) cycle ! skip the jloop

                  ! check if the element has be colored
                  if (Colors(Neighbors(j)) /= 0) then
                      ! add this color to blocked
                      BlockedColors(j) = Colors(Neighbors(j))
                  end if

              end do !j

              ! loop through existing number of colors
              do j = 1, NColor 
                  ! check if colorcount has reached nthreads
                  if (ColorCount(j) /= nthreads) then
                      if (allocated(tmp)) deallocate(tmp)
                      tmp  = findloc(BlockedColors, j)
                      isfree = size(tmp)==1 .and. tmp(1)==0
                      if (isfree) then
                          ! color elemment i with Color j
                          Colors(i) = j
                          ColorCount(j) = ColorCount(j) + 1
                          exit
                      end if
                  end if
              end do

              ! element i has not been colored
              ! either color(i) being blocked or have reached nthreads
              ! Color element i by a new color
              if (Colors(i)==0) then
                  NColor = NColor + 1
                  Colors(i) = NColor
                  ColorCount(NColor) = 1
              end if
          end do ! i
      end subroutine ColorByNeighbors 
END PROGRAM

