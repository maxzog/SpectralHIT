program hw4_q1
   use spectral
   use grid_class
   use navierstokes_class
   implicit none

   ! Variables
   integer :: N
   logical :: alias

   ! Command-line arguments
   integer :: iarg
   character(len=20) :: arg
   real(8) :: visc

   ! Grid type
   type(grid) :: mesh

   ! State type
   type(problem) :: state

   ! Size
   ! integer, parameter :: N=32
   
   ! State 
   double complex, dimension(:), allocatable :: u
   double complex, allocatable :: tempArray(:,:,:)

   ! Defaults
   N = 32                  ! Default grid size
   alias = .true.          ! Default to dealiasing enabled
   visc=1.0d0

   ! Read command-line arguments
   do iarg = 1, command_argument_count()
      call get_command_argument(iarg, arg)
      select case (iarg)
         case (1)
            read(arg, *) N
         case (2)
            read(arg, *) alias
         case (3)
            read(arg, *) visc
      end select
   end do

   ! Print the inputs for confirmation
   print *, 'Grid size N:', N
   print *, 'Alias:', alias
   print *, 'Visc:', visc

   allocate(u(N))

   ! Initialize domain
   mesh=grid(Nx=N,Ny=N,Nz=N,Lx=2*PI,Ly=2*PI,Lz=2*PI)
   call mesh%print

   state=problem_constructor(mesh=mesh, scheme=RK4, alias=alias, visc=visc)

   init_timer: block
      real(8) :: dt_max
      ! Calculate the maximum allowable dt for stability
      state%tf=0.5
      state%t=0.0
      state%dt=0.0001
      state%step=0
      state%stepf=99999999
   end block init_timer

   !> Init Taylor-Green vortex
   init_problem: block
      integer :: i,j,k
      do k=1,mesh%nz
         do j=1,mesh%ny
            do i=1,mesh%nx
               state%U(1,i,j,k)=-cos(mesh%xv(i))*sin(mesh%yv(j))
               state%U(2,i,j,k)=sin(mesh%xv(i))*cos(mesh%yv(j))
               state%U(3,i,j,k)=0.0d0
               state%P(i,j,j)=-0.25d0*(cos(2*mesh%xv(i)) + cos(2*mesh%yv(j)))
            end do
         end do
      end do
   end block init_problem

   call write_complex_array_3D(state%U(1,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/U_ic.txt")
   call write_complex_array_3D(state%U(2,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/V_ic.txt")
   call write_complex_array_3D(state%U(3,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/W_ic.txt")

   ! Transform the initial condition
   call state%transform_vel(mesh, FORWARD)

   ! call write_complex_array_3D(state%U, mesh%nx, "./outs/state_ic.txt")

   ! Advance Fourier coeffs
   do while(.not.(state%sdone.or.state%tdone))
      print *, "Step :: ", state%step
      call state%advance(mesh)
      state%U(3,:,:,:)=0.0d0
      call state%adjust_time()
   end do

   ! Transform solution back to physical space
   call state%transform_vel(mesh, BACKWARD)

   call write_complex_array_3D(state%U(1,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/U_001.txt")
   call write_complex_array_3D(state%U(2,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/V_001.txt")
   call write_complex_array_3D(state%U(3,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/W_001.txt")

   ! call output(state, mesh)

   contains

   subroutine output(state, mesh)
      implicit none
      type(problem), intent(in) :: state
      type(grid), intent(in) :: mesh
      integer :: i

      character(len=100) :: filename
      character(len=10) :: n_dir, visc_dir
      character(len=5) :: n_str, visc_str
   
      ! Convert N and visc to character strings
      write(n_str, '(I0)') N
      write(visc_str, '(F3.1)') visc
   
      ! Create directory names based on N and viscosity
      write(n_dir, '(A)') "n" // trim(adjustl(n_str))
      write(visc_dir, '(A)') "nu" // trim(adjustl(visc_str)) 
   
      ! Construct the filename based on alias setting
      if (state%dealias) then
         write(filename, '(A,A,A,A,A,I5.5,A)') "./", trim(n_dir), "/", trim(visc_dir), "/outs/state_rk4_", state%step, ".txt"
      else
         write(filename, '(A,A,A,A,A,I5.5,A)') "./", trim(n_dir), "/", trim(visc_dir), "/outs_no/state_rk4_", state%step, ".txt"
      end if

      call system("mkdir -p " // trim(n_dir) // "/" // trim(visc_dir) // "/outs")
      call system("mkdir -p " // trim(n_dir) // "/" // trim(visc_dir) // "/outs_no")
   
      ! Write the complex array to the generated filename
      call write_complex_array_1D(state%U, mesh%nx, filename)
   
      ! Print the current time and step
      print *, "Time :: ", state%t, " Step :: ", state%step
   end subroutine output

end program hw4_q1
