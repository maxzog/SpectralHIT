module grid_class
   implicit none

   public :: grid

   type grid
      integer :: dim
      integer :: Nx, Ny, Nz
      real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
      real(8) :: dx,dy,dz
      real(8) :: Lx,Ly,Lz

      real(8), dimension(:), allocatable :: xv,yv,zv
      integer, dimension(:), allocatable :: k1v,k2v,k3v
      
      contains

      procedure :: print => print_grid_info

      ! procedure :: solve_poisson
      procedure :: solve_poisson_1D
      procedure :: solve_poisson_2D
      procedure :: solve_poisson_3D

   end type grid

   interface grid
      procedure constructor
   end interface grid


   contains

   function constructor(Nx,Ny,Nz,Lx,Ly,Lz) result(self)
      implicit none
      type(grid) :: self
      real(8), optional :: Ly, Lz
      real(8) :: Lx
      integer, optional :: Ny, Nz
      integer :: Nx,i

      ! Start assuming 1D
      self%dim=1

      ! Set size of grid
      self%Nx=Nx
      if (present(Ny)) then
         self%Ny=Ny
         self%dim=self%dim+1
      else
         self%Ny=1
      end if
      if (present(Nz)) then
         self%Nz=Nz
         self%dim=self%dim+1
      else
         self%Nz=1
      end if

      ! Set physical dimensions
      self%Lx=Lx
      if (present(Ly)) then
         self%Ly=Ly
      else
         self%Ly=0.0
      end if
      if (present(Lz)) then
         self%Lz=Lz
      else
         self%Lz=0.0
      end if

      self%dx=self%Lx/self%Nx
      self%dy=self%Ly/self%Ny
      self%dz=self%Lz/self%Nz

      self%xmin=0.0; self%xmax=self%Lx
      self%ymin=0.0; self%ymax=self%Ly
      self%zmin=0.0; self%zmax=self%Lz

      allocate(self%xv(1:self%Nx+1)); self%xv=0.0
      if (self%dim.gt.1) allocate(self%yv(1:self%Ny+1))
      if (self%dim.gt.2) allocate(self%zv(1:self%Nz+1))

      allocate(self%k1v(1:self%Nx)); self%k1v=0 
      if (self%dim.gt.1) then
         allocate(self%k2v(1:self%Ny))
         self%k2v=0
      end if
      if (self%dim.gt.2) then
         allocate(self%k3v(1:self%Nz))
         self%k3v=0
      end if

      do i=1,self%Nx+1
         self%xv(i)=(i-1)*self%dx
      end do
      if (self%dim.gt.1) then
         do i=1,self%Ny+1
            self%yv(i)=(i-1)*self%dy
         end do 
      end if
      if (self%dim.gt.2) then
         do i=1,self%Nz+1
            self%zv(i)=(i-1)*self%dz
         end do
      end if
   end function constructor

   ! subroutine solve_poisson(this,k1v,k2v,k3v)
   !    implicit none
   !    class(grid), intent(inout) :: this
   !    integer, dimension(:), intent(in) :: k1v
   !    integer, dimension(:), intent(in), optional :: k2v, k3v

   !    if (this%dim.eq.1) call this%solve_poisson_1D(k1v)
   !    if (this%dim.eq.2) call this%solve_poisson_2D(k1v, k2v)
   !    if (this%dim.eq.3) call this%solve_poisson_3D(k1v, k2v, k3v)
   ! end subroutine solve_poisson

   subroutine solve_poisson_1D(this, k1v)
      implicit none
      class(grid), intent(in) :: this
      integer, dimension(:), intent(in) :: k1v
      ! TODO: Write this
   end subroutine solve_poisson_1D

   subroutine solve_poisson_2D(this, P, Q)
      implicit none
      class(grid), intent(inout) :: this
      complex(8), dimension(:,:), intent(inout) :: P, Q
      real(8) :: bottom
      integer :: i, j

      ! Solve
      do j=1,this%ny
         do i=1,this%nx
            bottom=sqrt(real(this%k1v(i)**2+this%k2v(j)**2,8))
            if (bottom.gt.0) then
               P(i,j)=-Q(i,j)/bottom
            else
               P(i,j)=0
            end if
         end do
      end do
   end subroutine solve_poisson_2D

   subroutine solve_poisson_3D(this, k1v, k2v, k3v)
      implicit none
      class(grid), intent(in) :: this
      integer, dimension(:), intent(in) :: k1v, k2v, k3v
      ! TODO: Write this
   end subroutine solve_poisson_3D

   subroutine print_grid_info(this)
     class(grid), intent(in) :: this
 
     ! Print header
     print *, "=============================="
     print *, "      Mesh Information         "
     print *, "=============================="
 
     ! Print dimension info
     print *, "Dimension: ", this%dim
     print *, ""
 
     ! Print number of points and domain ranges
     print *, "Number of points in each direction:"
     if (this%dim >= 1) print *, "  Nx =", this%Nx, "    Range: [", this%xmin, ",", this%xmax, "]"
     if (this%dim >= 2) print *, "  Ny =", this%Ny, "    Range: [", this%ymin, ",", this%ymax, "]"
     if (this%dim == 3) print *, "  Nz =", this%Nz, "    Range: [", this%zmin, ",", this%zmax, "]"
     print *, ""
 
     ! Print grid spacings
     print *, "Grid spacing in each direction:"
     if (this%dim >= 1) print *, "  dx =", this%dx
     if (this%dim >= 2) print *, "  dy =", this%dy
     if (this%dim == 3) print *, "  dz =", this%dz
     print *, ""
 
     ! Print domain lengths
     print *, "Domain lengths:"
     if (this%dim >= 1) print *, "  Lx =", this%Lx
     if (this%dim >= 2) print *, "  Ly =", this%Ly
     if (this%dim == 3) print *, "  Lz =", this%Lz
     print *, ""
 
     ! If grid points are allocated, print size of the arrays
     if (allocated(this%xv)) print *, "Grid points allocated in x-direction, size: ", size(this%xv)
     if (allocated(this%yv)) print *, "Grid points allocated in y-direction, size: ", size(this%yv)
     if (allocated(this%zv)) print *, "Grid points allocated in z-direction, size: ", size(this%zv)
 
     ! Print footer
     print *, "=============================="
     print *, ""
   end subroutine print_grid_info

end module grid_class
