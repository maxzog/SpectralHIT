module navierstokes_class
   use grid_class, only: grid

   public :: problem

   !> Integration scheme
   integer, parameter :: euler = 1
   integer, parameter :: RK2 = 2
   integer, parameter :: RK4 = 3

   !> Transform flags
   integer, parameter :: FORWARD = 1
   integer, parameter :: BACKWARD = 2

   type problem
      double complex, dimension(:, :, :, :), allocatable :: U
      double complex, dimension(:, :, :, :), allocatable :: U2
      !!      (uu uv uw)   (1 4 6)
      !! U2 = (uv vv vw) = (4 2 5)
      !!      (uw vw ww)   (6 5 3)
      double complex, dimension(:, :, :, :), allocatable :: RHS
      double complex, dimension(:, :, :), allocatable :: P
      real(8) :: visc = 0.0
      real(8) :: c = 0.0
      integer :: scheme = euler
      real(8) :: CFL = 0.0
      real(8) :: dt = 0.0
      real(8) :: tf = 0.0
      real(8) :: t = 0.0
      integer :: step = 0
      integer :: stepf = 0
      logical :: dealias = .true.
      logical :: sdone = .false.
      logical :: tdone = .false.
   contains
      procedure, public :: advance
      procedure, public :: get_rhs
      procedure, public :: transform_vel
      procedure, public :: transform_U2
      procedure, public :: compute_U2
      procedure, public :: pad_and_transform_U2
      procedure, public :: compute_phat
      procedure, public :: adjust_time
      procedure, public :: initialize_hit
   end type problem

   interface problem
      procedure problem_constructor
   end interface problem

contains

   function problem_constructor(mesh, scheme, alias, visc) result(self)
      implicit none
      type(problem) :: self
      type(grid) :: mesh
      integer :: scheme
      real(8), optional :: visc
      logical, optional :: alias

      allocate (self%P(1:mesh%nx, 1:mesh%ny, 1:mesh%nz)); self%P = 0.0
      allocate (self%U(3, 1:mesh%nx, 1:mesh%ny, 1:mesh%nz)); self%U = 0.0
      allocate (self%U2(6, 1:mesh%nx, 1:mesh%ny, 1:mesh%nz)); self%U2 = 0.0
      allocate (self%RHS(3, 1:mesh%nx, 1:mesh%ny, 1:mesh%nz)); self%RHS = 0.0

      if (present(alias)) self%dealias = alias

      self%scheme = scheme

      if (present(visc)) then
         self%visc = visc
      else
         self%visc = 0.0
      end if

      ! call random_init(.false., .false.)
      call random_init(.true., .true.)
   end function problem_constructor

   subroutine advance(this, mesh)
      implicit none
      class(problem), intent(inout) :: this
      class(grid), intent(in) :: mesh

      select case (this%scheme)
      case (euler)
         call advance_euler(this, mesh)
      case (RK2)
         call advance_rk2(this, mesh)
      case (RK4)
         call advance_rk4(this, mesh)
      end select

   contains

      subroutine advance_euler(this, mesh)
         implicit none
         class(problem), intent(inout) :: this
         class(grid), intent(in) :: mesh

         call this%get_rhs(mesh)
         this%U = this%U + this%RHS*this%dt
      end subroutine advance_euler

      subroutine advance_rk2(this, mesh)
         implicit none
         class(problem), intent(inout) :: this
         class(grid), intent(in) :: mesh
         double complex, dimension(:, :, :, :), allocatable :: Utemp

         allocate (Utemp, MOLD=this%U); Utemp = this%U

         call this%get_rhs(mesh)
         this%U = this%U + 0.5d0*this%RHS*this%dt
         call this%get_rhs(mesh)
         this%U = Utemp + this%RHS*this%dt

         deallocate (Utemp)
      end subroutine advance_rk2

      subroutine advance_rk4(this, mesh)
         implicit none
         class(problem), intent(inout) :: this
         class(grid), intent(in) :: mesh
         double complex, dimension(:, :, :, :), allocatable :: Utemp
         double complex, dimension(:, :, :, :), allocatable :: k1, k2, k3, k4

         ! Store state at step n
         allocate (Utemp, MOLD=this%U); Utemp = this%U

         ! Create storage for intermediate steps
         allocate (k1, MOLD=this%U); k1 = 0.0d0
         allocate (k2, MOLD=this%U); k2 = 0.0d0
         allocate (k3, MOLD=this%U); k3 = 0.0d0

         ! Call RHS for first stage
         call this%get_rhs(mesh)
         k1 = this%RHS
         ! Advance using first stage
         this%U = Utemp + 0.5d0*this%dt*k1
         call this%get_rhs(mesh)
         k2 = this%RHS
         ! Advance using second stage
         this%U = Utemp + 0.5d0*this%dt*k2
         call this%get_rhs(mesh)
         k3 = this%RHS
         this%U = Utemp + this%dt*k3
         ! Call RHS for fourth stage
         call this%get_rhs(mesh)

         this%U = Utemp + this%dt/6.0d0*(k1 + 2.0d0*k2 + 2.0d0*k3 + this%RHS)

         deallocate (Utemp, k1, k2, k3)
      end subroutine advance_rk4
   end subroutine advance

   subroutine get_rhs(this, mesh)
      implicit none
      class(problem), intent(inout) :: this
      class(grid), intent(in) :: mesh

      call get_rhs_ns(this)

   contains
      !> 3D Navier-Stokes RHS without dealiasing
      subroutine get_rhs_ns(this)
         use spectral
         implicit none
         class(problem), intent(inout) :: this
         double complex :: duu
         integer :: i, j, k
         real(8) :: kmag

         if (this%dealias) then
            ! ! Compute non-linear term
            ! call this%compute_U2(mesh)
            ! Pad and transform U2
            call this%pad_and_transform_U2(mesh)
         else
            ! Bring velocity to physical space
            call this%transform_vel(mesh, BACKWARD)
            ! Compute non-linear term
            call this%compute_U2(mesh)
            ! Bring non-linear term to spectral space
            call this%transform_U2(mesh, FORWARD)
            ! Bring velocity back to spectral space
            call this%transform_vel(mesh, FORWARD)
         end if
         ! Compute pressure in spectral space
         call this%compute_phat(mesh)

         !$omp parallel do private(i, j, k, kmag)
         do k = 1, mesh%Nz
            do j = 1, mesh%Ny
               do i = 1, mesh%Nx
                  kmag = sqrt(real(mesh%k1v(i)**2 + mesh%k2v(j)**2 + mesh%k3v(k)**2, 8))
                  ! RHSi = - d/dxj (uiuj) + d/dxi P + nu d^2/dxj^2 u_i
                  ! RHS - x
                  duu = i_unit*(mesh%k1v(i)*this%U2(1, i, j, k) + mesh%k2v(j)*this%U2(4, i, j, k) + mesh%k3v(k)*this%U2(6, i, j, k))
                  this%RHS(1, i, j, k) = -duu - i_unit*mesh%k1v(i)*this%P(i, j, k) - this%visc*kmag**2*this%U(1, i, j, k)
                  ! RHS - y
                  duu = i_unit*(mesh%k1v(i)*this%U2(4, i, j, k) + mesh%k2v(j)*this%U2(2, i, j, k) + mesh%k3v(k)*this%U2(5, i, j, k))
                  this%RHS(2, i, j, k) = -duu - i_unit*mesh%k2v(j)*this%P(i, j, k) - this%visc*kmag**2*this%U(2, i, j, k)
                  ! RHS - z
                  duu = i_unit*(mesh%k1v(i)*this%U2(6, i, j, k) + mesh%k2v(j)*this%U2(5, i, j, k) + mesh%k3v(k)*this%U2(3, i, j, k))
                  this%RHS(3, i, j, k) = -duu - i_unit*mesh%k3v(k)*this%P(i, j, k) - this%visc*kmag**2*this%U(3, i, j, k)
               end do
            end do
         end do
         !$omp end parallel do
      end subroutine get_rhs_ns
   end subroutine get_rhs

   !> Compute the pressure in spectral space
   subroutine compute_phat(this, mesh)
      implicit none
      class(problem), intent(inout) :: this
      class(grid), intent(in) :: mesh
      integer :: i, j, k
      real(8) :: kmag

      this%P = 0.0
      !$omp parallel do private(i, j, k, kmag)
      do k = 1, mesh%nz
         do j = 1, mesh%ny
            do i = 1, mesh%nx
               kmag = sqrt(real(mesh%k1v(i)**2 + mesh%k2v(j)**2 + mesh%k3v(k)**2, 8))
               if (kmag .lt. epsilon(1.0d0)) cycle
               this%P(i, j, k) = this%P(i, j, k) - mesh%k1v(i)**2*this%U2(1, i, j, k)
               this%P(i, j, k) = this%P(i, j, k) - mesh%k2v(j)**2*this%U2(2, i, j, k)
               this%P(i, j, k) = this%P(i, j, k) - mesh%k3v(k)**2*this%U2(3, i, j, k)
               this%P(i, j, k) = this%P(i, j, k) - 2.0*mesh%k1v(i)*mesh%k2v(j)*this%U2(4, i, j, k)
               this%P(i, j, k) = this%P(i, j, k) - 2.0*mesh%k2v(j)*mesh%k3v(k)*this%U2(5, i, j, k)
               this%P(i, j, k) = this%P(i, j, k) - 2.0*mesh%k1v(i)*mesh%k3v(k)*this%U2(6, i, j, k)
               this%P(i, j, k) = this%P(i, j, k)/kmag**2
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine compute_phat

   !> Compute stress tensor
   subroutine compute_U2(this, mesh)
      implicit none
      class(problem), intent(inout) :: this
      class(grid), intent(in) :: mesh
      integer :: i, j, k

      ! Compute non-linear term in physical space
      !$omp parallel do private(i, j, k)
      do k = 1, mesh%Nz
         do j = 1, mesh%Ny
            do i = 1, mesh%Nx
               this%U2(1, i, j, k) = this%U(1, i, j, k)*this%U(1, i, j, k)
               this%U2(2, i, j, k) = this%U(2, i, j, k)*this%U(2, i, j, k)
               this%U2(3, i, j, k) = this%U(3, i, j, k)*this%U(3, i, j, k)

               this%U2(4, i, j, k) = this%U(1, i, j, k)*this%U(2, i, j, k)
               this%U2(5, i, j, k) = this%U(2, i, j, k)*this%U(3, i, j, k)
               this%U2(6, i, j, k) = this%U(1, i, j, k)*this%U(3, i, j, k)
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine compute_U2

   !> Transform U2 with padding to prevent aliasing
   subroutine pad_and_transform_U2(this, mesh)
      use spectral
      implicit none
      class(problem), intent(inout) :: this
      class(grid), intent(in) :: mesh
      double complex, dimension(:, :, :, :), allocatable :: U2pad, Upad
      double complex, dimension(:, :, :), allocatable   :: tempArr
      integer :: i, j, k
      logical, allocatable :: isMiddleX(:), isMiddleY(:), isMiddleZ(:)
      integer, allocatable :: shiftX(:), shiftY(:), shiftZ(:)

      allocate (Upad(3, 1:2*mesh%nx, 1:2*mesh%ny, 1:2*mesh%nz))
      Upad = cmplx(0.0d0, 0.0d0)

      allocate (U2pad(6, 1:2*mesh%nx, 1:2*mesh%ny, 1:2*mesh%nz))
      U2pad = cmplx(0.0d0, 0.0d0)

      allocate (tempArr(1:2*mesh%nx, 1:2*mesh%ny, 1:2*mesh%nz))
      tempArr = cmplx(0.0d0, 0.0d0)

      allocate (shiftX(1:2*mesh%nx), shiftY(1:2*mesh%ny), shiftZ(1:2*mesh%nz))
      allocate (isMiddleX(1:2*mesh%nx), isMiddleY(1:2*mesh%ny), isMiddleZ(1:2*mesh%nz))

      ! Precompute shifts and middle-region masks
      do i = 1, 2*mesh%nx
         shiftX(i) = merge(1, 0, i > 1.5*mesh%nx)
         isMiddleX(i) = i > 0.5*mesh%nx .and. i <= 1.5*mesh%nx
      end do

      do j = 1, 2*mesh%ny
         shiftY(j) = merge(1, 0, j > 1.5*mesh%ny)
         isMiddleY(j) = j > 0.5*mesh%ny .and. j <= 1.5*mesh%ny
      end do

      do k = 1, 2*mesh%nz
         shiftZ(k) = merge(1, 0, k > 1.5*mesh%nz)
         isMiddleZ(k) = k > 0.5*mesh%nz .and. k <= 1.5*mesh%nz
      end do

      ! Fill U2pad
      !$omp parallel do private(i, j, k)
      do k = 1, 2*mesh%nz
         do j = 1, 2*mesh%ny
            do i = 1, 2*mesh%nx
               if (.not. (isMiddleX(i) .or. isMiddleY(j) .or. isMiddleZ(k))) then
                  Upad(:, i, j, k) = this%U(:, i - shiftX(i)*mesh%nx, j - shiftY(j)*mesh%ny, k - shiftZ(k)*mesh%nz)
               end if
            end do
         end do
      end do
      !$omp end parallel do

      ! Transform the padded U array
      do i = 1, 3
         tempArr = Upad(i, :, :, :)
         call iFFT_3D(tempArr, tempArr, 2*mesh%nx, 2*mesh%ny, 2*mesh%nz)
         Upad(i, :, :, :) = tempArr
      end do

      ! Compute U2 in physical space
      do k=1,2*mesh%nz
         do j=1,2*mesh%ny
            do i=1,2*mesh%nx
               U2pad(1, i, j, k) = Upad(1, i, j, k)*Upad(1, i, j, k)
               U2pad(2, i, j, k) = Upad(2, i, j, k)*Upad(2, i, j, k)
               U2pad(3, i, j, k) = Upad(3, i, j, k)*Upad(3, i, j, k)

               U2pad(4, i, j, k) = Upad(1, i, j, k)*Upad(2, i, j, k)
               U2pad(5, i, j, k) = Upad(2, i, j, k)*Upad(3, i, j, k)
               U2pad(6, i, j, k) = Upad(1, i, j, k)*Upad(3, i, j, k)
            end do
         end do
      end do

      ! Transform the padded U2 array
      do i = 1, 6
         tempArr = U2pad(i, :, :, :)
         call FFT_3D(tempArr, tempArr, 2*mesh%nx, 2*mesh%ny, 2*mesh%nz)
         U2pad(i, :, :, :) = tempArr
      end do

      ! Copy back transformed data
      !$omp parallel do private(i, j, k)
      do k = 1, 2*mesh%nz
         do j = 1, 2*mesh%ny
            do i = 1, 2*mesh%nx
               if (.not. (isMiddleX(i) .or. isMiddleY(j) .or. isMiddleZ(k))) then
                  this%U2(:, i - shiftX(i)*mesh%nx, j - shiftY(j)*mesh%ny, k - shiftZ(k)*mesh%nz) = U2pad(:, i, j, k)
               end if
            end do
         end do
      end do
      !$omp end parallel do

      deallocate (U2pad, Upad, tempArr, shiftX, shiftY, shiftZ, isMiddleX, isMiddleY, isMiddleZ)
   end subroutine pad_and_transform_U2

   ! !> Transform U2 with padding to prevent aliasing
   ! subroutine pad_and_transform_U2(this, mesh)
   !    use spectral
   !    implicit none
   !    class(problem), intent(inout) :: this
   !    class(grid), intent(in) :: mesh
   !    double complex, dimension(:, :, :, :), allocatable :: U2pad
   !    double complex, dimension(:, :, :), allocatable   :: tempArr
   !    integer :: i, j, k
   !    logical, allocatable :: isMiddleX(:), isMiddleY(:), isMiddleZ(:)
   !    integer, allocatable :: shiftX(:), shiftY(:), shiftZ(:)

   !    allocate (U2pad(6, 1:2*mesh%nx, 1:2*mesh%ny, 1:2*mesh%nz))
   !    U2pad = cmplx(0.0d0, 0.0d0)

   !    allocate (tempArr(1:2*mesh%nx, 1:2*mesh%ny, 1:2*mesh%nz))
   !    tempArr = cmplx(0.0d0, 0.0d0)

   !    allocate (shiftX(1:2*mesh%nx), shiftY(1:2*mesh%ny), shiftZ(1:2*mesh%nz))
   !    allocate (isMiddleX(1:2*mesh%nx), isMiddleY(1:2*mesh%ny), isMiddleZ(1:2*mesh%nz))

   !    ! Precompute shifts and middle-region masks
   !    do i = 1, 2*mesh%nx
   !       shiftX(i) = merge(1, 0, i > 1.5*mesh%nx)
   !       isMiddleX(i) = i > 0.5*mesh%nx .and. i <= 1.5*mesh%nx
   !    end do

   !    do j = 1, 2*mesh%ny
   !       shiftY(j) = merge(1, 0, j > 1.5*mesh%ny)
   !       isMiddleY(j) = j > 0.5*mesh%ny .and. j <= 1.5*mesh%ny
   !    end do

   !    do k = 1, 2*mesh%nz
   !       shiftZ(k) = merge(1, 0, k > 1.5*mesh%nz)
   !       isMiddleZ(k) = k > 0.5*mesh%nz .and. k <= 1.5*mesh%nz
   !    end do

   !    ! Fill U2pad
   !    !$omp parallel do private(i, j, k)
   !    do k = 1, 2*mesh%nz
   !       do j = 1, 2*mesh%ny
   !          do i = 1, 2*mesh%nx
   !             if (.not. (isMiddleX(i) .or. isMiddleY(j) .or. isMiddleZ(k))) then
   !                U2pad(:, i, j, k) = this%U2(:, i - shiftX(i)*mesh%nx, j - shiftY(j)*mesh%ny, k - shiftZ(k)*mesh%nz)
   !             end if
   !          end do
   !       end do
   !    end do
   !    !$omp end parallel do

   !    call write_complex_array_3D(this%U2(1,:,:,:), mesh%nx, mesh%ny, mesh%nz, "./outs/U2.txt")
   !    call write_complex_array_3D(U2pad(1,:,:,:), 2*mesh%nx, 2*mesh%ny, 2*mesh%nz, "./outs/U2pad.txt")
   !    call abort

   !    ! Transform the padded U2 array
   !    do i = 1, 6
   !       tempArr = U2pad(i, :, :, :)
   !       call FFT_3D(tempArr, tempArr, 2*mesh%nx, 2*mesh%ny, 2*mesh%nz)
   !       U2pad(i, :, :, :) = tempArr
   !    end do

   !    ! Copy back transformed data
   !    !$omp parallel do private(i, j, k)
   !    do k = 1, 2*mesh%nz
   !       do j = 1, 2*mesh%ny
   !          do i = 1, 2*mesh%nx
   !             if (.not. (isMiddleX(i) .or. isMiddleY(j) .or. isMiddleZ(k))) then
   !                this%U2(:, i - shiftX(i)*mesh%nx, j - shiftY(j)*mesh%ny, k - shiftZ(k)*mesh%nz) = U2pad(:, i, j, k)
   !             end if
   !          end do
   !       end do
   !    end do
   !    !$omp end parallel do

   !    deallocate (U2pad, tempArr, shiftX, shiftY, shiftZ, isMiddleX, isMiddleY, isMiddleZ)
   ! end subroutine pad_and_transform_U2

   !> For convenience
   !> Performs three 3D transforms on the velocity field
   subroutine transform_vel(this, mesh, direction)
      use spectral
      implicit none
      class(problem), intent(inout) :: this
      class(grid), intent(in) :: mesh
      double complex, dimension(:, :, :), allocatable :: tempArr
      integer :: i, direction

      allocate (tempArr, MOLD=this%U(1, :, :, :))

      select case (direction)
      case (FORWARD)
         do i = 1, 3
            tempArr = this%U(i, :, :, :)
            call FFT_3D(tempArr, tempArr, mesh%nx, mesh%ny, mesh%nz)
            this%U(i, :, :, :) = tempArr
         end do
      case (BACKWARD)
         do i = 1, 3
            tempArr = this%U(i, :, :, :)
            call iFFT_3D(tempArr, tempArr, mesh%nx, mesh%ny, mesh%nz)
            this%U(i, :, :, :) = tempArr
         end do
      end select

      deallocate (tempArr)
   end subroutine transform_vel

   !> For convenience
   !> Performs six 3D transforms on the stress tensor
   subroutine transform_U2(this, mesh, direction)
      use spectral
      implicit none
      class(problem), intent(inout) :: this
      class(grid), intent(in) :: mesh
      double complex, dimension(:, :, :), allocatable :: tempArr
      integer :: i, direction

      allocate (tempArr, MOLD=this%U2(1, :, :, :))
      select case (direction)
      case (FORWARD)
         do i = 1, 6
            tempArr = this%U2(i, :, :, :)
            call FFT_3D(tempArr, tempArr, mesh%nx, mesh%ny, mesh%nz)
            this%U2(i, :, :, :) = tempArr
         end do
      case (BACKWARD)
         do i = 1, 6
            tempArr = this%U2(i, :, :, :)
            call iFFT_3D(tempArr, tempArr, mesh%nx, mesh%ny, mesh%nz)
            this%U2(i, :, :, :) = tempArr
         end do
      end select
      deallocate (tempArr)
   end subroutine transform_U2

   subroutine initialize_hit(this, mesh, k0, u0)
      implicit none
      class(problem), intent(inout) :: this
      class(grid), intent(in) :: mesh
      real(8), intent(in) :: k0, u0
      double complex :: top, bot, alpha, beta
      integer :: i, j, k
      integer, dimension(3) :: kv
      real(8) :: kmag

      this%U = cmplx(0.0d0, 0.0d0)

      ! Set zeroth mode and oddball values
      this%U(:, 1, 1, 1) = cmplx(0.0d0, 0.0d0)
      do k = 1, mesh%Nz
         do j = 1, mesh%Ny
            do i = 1, mesh%nx
               if (i .eq. mesh%Nx/2 + 1 .or. j .eq. mesh%Ny/2 + 1 .or. k .eq. mesh%Nz/2 + 1) then
                  this%U(:, i, j, k) = cmplx(0.0d0, 0.0d0)  ! Ensure these are real
                  cycle
               end if
            end do
         end do
      end do

      do k = 2, mesh%nz
         do j = 2, mesh%ny
            do i = 2, mesh%nx/2
               ! Compute magnitude of wavenumber vector
               kmag = sqrt(real(mesh%k1v(i)**2 + mesh%k2v(j)**2 + mesh%k3v(k)**2, 8))

               if (kmag .lt. 10.0d0*epsilon(1.0d0)) cycle

               ! Get alpha and beta coefficients
               call get_alpha_beta(kmag, k0, u0, alpha, beta)

               ! u (x) component
               top = alpha*kmag*mesh%k2v(j) + beta*mesh%k1v(i)*mesh%k3v(k)
               bot = kmag*sqrt(real(mesh%k1v(i)**2 + mesh%k2v(j)**2, 8))
               if (abs(bot) > 10.0d0*epsilon(1.0d0)) then
                  this%U(1, i, j, k) = top/bot
                  this%U(1, mesh%nx - i + 2, mesh%ny - j + 2, mesh%nz - k + 2) = conjg(top/bot)
               else
                  this%U(1, i, j, k) = alpha
                  this%U(1, mesh%nx - i + 2, mesh%ny - j + 2, mesh%nz - k + 2) = conjg(alpha)
               end if

               ! v (y) component
               top = beta*mesh%k2v(j)*mesh%k3v(k) - alpha*kmag*mesh%k1v(i)
               bot = kmag*sqrt(real(mesh%k1v(i)**2 + mesh%k2v(j)**2, 8))
               if (abs(bot) .gt. 10.0d0*epsilon(1.0d0)) then
                  this%U(2, i, j, k) = top/bot
                  this%U(2, mesh%nx - i + 2, mesh%ny - j + 2, mesh%nz - k + 2) = conjg(top/bot)
               else
                  this%U(2, i, j, k) = beta
                  this%U(2, mesh%nx - i + 2, mesh%ny - j + 2, mesh%nz - k + 2) = conjg(beta)
               end if

               ! w (z) component
               top = -beta*sqrt(real(mesh%k1v(i)**2 + mesh%k2v(j)**2, 8))
               bot = kmag
               if (abs(top) .gt. 10.0d0*epsilon(1.0d0)) then
                  this%U(3, i, j, k) = top/bot
                  this%U(3, mesh%nx - i + 2, mesh%ny - j + 2, mesh%nz - k + 2) = conjg(top/bot)
               else
                  this%U(3, i, j, k) = cmplx(0.0d0, 0.0d0)
                  this%U(3, mesh%nx - i + 2, mesh%ny - j + 2, mesh%nz - k + 2) = cmplx(0.0d0, 0.0d0)
               end if
            end do
         end do
      end do
   end subroutine initialize_hit

   subroutine get_alpha_beta(kmag, k0, u0, alpha, beta)
      use spectral
      real(8), intent(in) :: kmag, k0, u0
      double complex, intent(inout) :: alpha, beta
      real(8) :: theta1, theta2, phi

      ! Get random angles
      call random_number(theta1)
      call random_number(theta2)
      call random_number(phi)

      ! Scale angles to be in [0,2pi) or [-pi,pi)
      theta1 = theta1*2.0*PI
      theta2 = theta2*2.0*PI
      theta1 = (theta1*0.5 - 1.0)*2.0*PI
      theta2 = (theta2*0.5 - 1.0)*2.0*PI
      phi = phi*2.0*PI

      ! Compute coeffs
      alpha = sqrt(spectrum(kmag, k0, u0)/(4*PI*kmag**2))*exp(i_unit*theta1)*cos(phi)
      beta = sqrt(spectrum(kmag, k0, u0)/(4*PI*kmag**2))*exp(i_unit*theta2)*sin(phi)
   end subroutine get_alpha_beta

   ! Prescribed energy spectrum for initial condition
   function spectrum(kmag, k0, u0) result(energy)
      use spectral
      real(8) :: energy, kmag, k0, u0
      energy = 16.0*sqrt(2.0/PI)*u0**2/k0*(kmag/k0)**4*exp(-2.0*(kmag/k0)**2)
   end function spectrum

   subroutine adjust_time(this)
      implicit none
      class(problem), intent(inout) :: this

      this%t = this%t + this%dt
      this%step = this%step + 1

      this%tdone = this%t > this%tf
      this%sdone = this%step > this%stepf
   end subroutine adjust_time

end module navierstokes_class
