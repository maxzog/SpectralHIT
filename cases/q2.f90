program hw5
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
   real(8) :: visc, k0, u0

   ! Type inits for grid and state 
   type(grid) :: mesh
   type(problem) :: state

   ! Energy spectrum
   real(8), dimension(:), allocatable :: Ek,Ekv

   ! Stats
   real(8) :: divergence, skewness, kinetic_energy
   integer :: unit_output
   character(len=100) :: output_file

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

   ! Prep spectrum and axis
   allocate(Ek (1:N)); Ek=0.0d0
   allocate(Ekv(1:N)); Ekv=0.0d0

   ! Initialize domain
   mesh=grid(Nx=N,Ny=N,Nz=N,Lx=2*PI,Ly=2*PI,Lz=2*PI)
   mesh%k1v=get_kv(mesh%Nx)
   mesh%k2v=get_kv(mesh%Ny)
   mesh%k3v=get_kv(mesh%Nz)

   state=problem_constructor(mesh=mesh, scheme=RK4, alias=alias, visc=visc)

   init_timer: block
      real(8) :: dt_max
      ! Calculate the maximum allowable dt for stability
      state%tf=5.0
      state%t=0.0
      state%dt=0.001
      state%step=0
      state%stepf= 99999999
   end block init_timer

   k0=5.0d0
   u0=1.0d0

   call state%initialize_hit(mesh,k0,u0)

   state%U=state%U*PI**(0.25d0)

   call write_complex_array_3D(state%U(1,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/Uh_ic.txt")
   call write_complex_array_3D(state%U(2,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/Vh_ic.txt")
   call write_complex_array_3D(state%U(3,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/Wh_ic.txt")

   call compute_radial_spectrum(state%U, mesh, Ek, Ekv)

   call write_double_array(Ek, N, "./outs/Ek_ic.txt")
   call write_double_array(Ekv, N, "./outs/Ekax_ic.txt")

   ! Transform the initial condition
   call state%transform_vel(mesh, BACKWARD)

   call write_complex_array_3D(state%U(1,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/U_ic.txt")
   call write_complex_array_3D(state%U(2,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/V_ic.txt")
   call write_complex_array_3D(state%U(3,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/W_ic.txt")

   ! Transform the initial condition
   call state%transform_vel(mesh, FORWARD)

   ! Advance Fourier coeffs
   do while(.not.(state%sdone.or.state%tdone))
      call state%advance(mesh)
      call state%adjust_time()
      ! Inside the loop
      if (mod(state%step,10) .eq. 1) then
         call compute_radial_spectrum(state%U, mesh, Ek, Ekv)
         call compute_stats(state, mesh, divergence, skewness)
         call state%transform_vel(mesh, BACKWARD)
         call compute_tke(state, mesh, kinetic_energy)
         call output(state, mesh)
         call state%transform_vel(mesh, FORWARD)
         output_file = 'output.txt'
         ! Open file
         open(newunit=unit_output, file=output_file, status='unknown', action='write', position='append')
         ! Write to the file
         write(unit_output, '(I10, 3F15.6)') state%step, divergence, skewness, kinetic_energy
         ! Close file
         close(unit_output)
      end if
   end do


   ! Transform the initial condition
   call state%transform_vel(mesh, BACKWARD)

   call write_complex_array_3D(state%U(1,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/U_later.txt")
   call write_complex_array_3D(state%U(2,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/V_later.txt")
   call write_complex_array_3D(state%U(3,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/W_later.txt")

   deallocate(Ek, Ekv)

   contains

   subroutine compute_tke(state, mesh, tke)
      implicit none
      type(problem), intent(in) :: state
      type(grid), intent(in) :: mesh
      real(8), intent(inout) :: tke 
      real(8) :: mx,my,mz,vx,vy,vz
      integer :: i,j,k

      ! Get mean velocity components
      mx=0.0d0
      my=0.0d0
      mz=0.0d0
      do k=1,mesh%Nz
         do j=1,mesh%Ny
            do i=1,mesh%Nx
               mx = mx + state%U(1,i,j,k)
               my = my + state%U(2,i,j,k)
               mz = mz + state%U(3,i,j,k)
            end do
         end do
      end do
      mx=mx/(mesh%nx*mesh%ny*mesh%nz)
      my=my/(mesh%nx*mesh%ny*mesh%nz)
      mz=mz/(mesh%nx*mesh%ny*mesh%nz)

      ! Compute variance
      vx=0.0d0
      vy=0.0d0
      vz=0.0d0
      do k=1,mesh%Nz
         do j=1,mesh%Ny
            do i=1,mesh%Nx
               vx = vx + (state%U(1,i,j,k) - mx)**2
               vy = vy + (state%U(2,i,j,k) - my)**2
               vz = vz + (state%U(3,i,j,k) - mz)**2
            end do
         end do
      end do
      vx=vx/(mesh%nx*mesh%ny*mesh%nz)
      vy=vy/(mesh%nx*mesh%ny*mesh%nz)
      vz=vz/(mesh%nx*mesh%ny*mesh%nz)

      ! Get the tke
      tke=0.5d0*(vx+vy+vz)

   end subroutine compute_tke

   subroutine compute_stats(state, mesh, meanU, skew)
      implicit none
      type(problem), intent(inout) :: state
      type(grid), intent(in) :: mesh
      real(8), intent(inout) :: meanU, skew
      real(8) :: meanU2, meanU3
      double complex, dimension(:,:,:), allocatable :: tempArr
      integer :: i, j, k

      ! Compute the terms du_i / dx_i in Fourier space and store in RHS array
      do k=1,mesh%Nz
         do j=1,mesh%Ny
            do i=1,mesh%Nx
               state%RHS(1,i,j,k) = i_unit*mesh%k1v(i)*state%U(1,i,j,k)
               state%RHS(2,i,j,k) = i_unit*mesh%k2v(j)*state%U(2,i,j,k)
               state%RHS(3,i,j,k) = i_unit*mesh%k3v(k)*state%U(3,i,j,k)
            end do
         end do
      end do

      allocate(tempArr, MOLD=state%RHS(1,:,:,:))

      ! Transform du_i / dx_i to physical space
      do i=1,3
         tempArr = state%RHS(i,:,:,:)
         call iFFT_3D(tempArr,tempArr,mesh%nx,mesh%ny,mesh%nz)
         state%RHS(i,:,:,:) = tempArr
      end do

      deallocate(tempArr)

      meanU =0.0d0
      meanU2=0.0d0
      meanU3=0.0d0
      do k=1,mesh%Nz
         do j=1,mesh%Ny
            do i=1,mesh%Nx
               meanU = meanU + sum(state%RHS(:,i,j,k))/3
               meanU2 = meanU2 + sum(state%RHS(:,i,j,k)**2)/3
               meanU3 = meanU3 + sum(state%RHS(:,i,j,k)**3)/3
            end do
         end do
      end do

      meanU=meanU/(mesh%nx*mesh%ny*mesh%nz)
      meanU2=meanU2/(mesh%nx*mesh%ny*mesh%nz)
      meanU3=meanU3/(mesh%nx*mesh%ny*mesh%nz)

      skew=meanU3/meanU2**(1.5d0)

   end subroutine compute_stats

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
   
      write(filename, '(A,I5.5,A)') "./outs/U_", state%step, ".txt"
      ! Write the complex array to the generated filename
      call write_complex_array_3D(state%U(1,:,:,:), mesh%nx, mesh%ny, mesh%nz, filename)
      write(filename, '(A,I5.5,A)') "./outs/V_", state%step, ".txt"
      ! Write the complex array to the generated filename
      call write_complex_array_3D(state%U(2,:,:,:), mesh%nx, mesh%ny, mesh%nz, filename)
      write(filename, '(A,I5.5,A)') "./outs/W_", state%step, ".txt"
      ! Write the complex array to the generated filename
      call write_complex_array_3D(state%U(3,:,:,:), mesh%nx, mesh%ny, mesh%nz, filename)

      write(filename, '(A,I5.5,A)') "./outs/Ek_", state%step, ".txt"
      call write_double_array(Ek, N, filename)
      write(filename, '(A,I5.5,A)') "./outs/Ekax_", state%step, ".txt"
      call write_double_array(Ekv, N, filename)
   
      ! Print the current time and step
      print *, "Time :: ", state%t, " Step :: ", state%step
   end subroutine output
end program hw5
