! orbfit_propagator_simple.f90
! Simplified OrbFit propagator for comparison with AstDyn
! Reads .eq1 file, propagates orbit, outputs Cartesian states

PROGRAM orbfit_propagator_simple
  IMPLICIT NONE
  
  ! Constants
  DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846d0
  DOUBLE PRECISION, PARAMETER :: DEG2RAD = PI / 180.d0
  DOUBLE PRECISION, PARAMETER :: AU = 1.49597870691d8  ! km
  DOUBLE PRECISION, PARAMETER :: GMS = 2.959122082855911d-4  ! AU^3/day^2
  
  ! Orbital elements (Keplerian)
  DOUBLE PRECISION :: a, e, i, Omega, omega, M, epoch_mjd
  CHARACTER(LEN=100) :: obj_name
  
  ! Cartesian state
  DOUBLE PRECISION, DIMENSION(6) :: state
  DOUBLE PRECISION :: t_start, t_end, dt, t_current
  INTEGER :: n_steps, step
  
  ! File I/O
  CHARACTER(LEN=256) :: eq1_file, output_file
  INTEGER :: unit_in, unit_out, ios
  
  ! Get command line arguments
  IF(COMMAND_ARGUMENT_COUNT() < 1) THEN
    WRITE(*,*) 'Usage: orbfit_propagator_simple <eq1_file> [output_file]'
    WRITE(*,*) 'Example: ./orbfit_propagator_simple 203_astdys.eq1 orbfit_states.csv'
    STOP
  END IF
  
  CALL GET_COMMAND_ARGUMENT(1, eq1_file)
  IF(COMMAND_ARGUMENT_COUNT() >= 2) THEN
    CALL GET_COMMAND_ARGUMENT(2, output_file)
  ELSE
    output_file = 'orbfit_states.csv'
  END IF
  
  ! Read orbital elements from .eq1 file
  CALL read_eq1_simple(TRIM(eq1_file), obj_name, epoch_mjd, a, e, i, Omega, omega, M)
  
  WRITE(*,'(A)') '=========================================='
  WRITE(*,'(A)') '  OrbFit Propagator (Simplified)'
  WRITE(*,'(A)') '=========================================='
  WRITE(*,'(A,A)') 'Object: ', TRIM(obj_name)
  WRITE(*,'(A,F12.4)') 'Epoch (MJD TDB): ', epoch_mjd
  WRITE(*,'(A,F12.6)') 'a (AU): ', a
  WRITE(*,'(A,F12.6)') 'e: ', e
  WRITE(*,'(A,F12.6)') 'i (deg): ', i / DEG2RAD
  
  ! Convert Keplerian to Cartesian at epoch
  CALL keplerian_to_cartesian(a, e, i, Omega, omega, M, GMS, state)
  
  ! Propagation parameters
  t_start = 61000.0d0  ! MJD
  t_end = 61030.0d0    ! MJD
  dt = 1.0d0           ! 1 day
  n_steps = NINT((t_end - t_start) / dt) + 1
  
  ! Open output file
  OPEN(NEWUNIT=unit_out, FILE=TRIM(output_file), STATUS='REPLACE', ACTION='WRITE')
  WRITE(unit_out,'(A)') 'MJD,X_AU,Y_AU,Z_AU,VX_AU_d,VY_AU_d,VZ_AU_d'
  
  ! Propagate and output
  WRITE(*,'(A,I0,A)') 'Propagating ', n_steps, ' steps...'
  
  DO step = 0, n_steps - 1
    t_current = t_start + step * dt
    
    ! Simple 2-body propagation (placeholder for full OrbFit propagator)
    ! In real implementation, this would call OrbFit's force model + integrator
    CALL propagate_twobody(state, epoch_mjd, t_current, GMS)
    
    ! Write state
    WRITE(unit_out,'(F12.4,6(",",ES23.15))') t_current, state(1:6)
  END DO
  
  CLOSE(unit_out)
  
  WRITE(*,'(A,A)') 'Output written to: ', TRIM(output_file)
  WRITE(*,'(A)') '=========================================='
  
CONTAINS

  ! Read simplified .eq1 file (Keplerian elements only)
  SUBROUTINE read_eq1_simple(filename, name, epoch, a, e, i, Omega, omega, M)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(OUT) :: name
    DOUBLE PRECISION, INTENT(OUT) :: epoch, a, e, i, Omega, omega, M
    
    INTEGER :: unit, ios
    CHARACTER(LEN=256) :: line
    LOGICAL :: in_data
    
    OPEN(NEWUNIT=unit, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF(ios /= 0) THEN
      WRITE(*,*) 'ERROR: Cannot open file: ', TRIM(filename)
      STOP
    END IF
    
    in_data = .FALSE.
    name = ''
    
    DO
      READ(unit,'(A)',IOSTAT=ios) line
      IF(ios /= 0) EXIT
      
      ! Skip header
      IF(INDEX(line, 'END_OF_HEADER') > 0) THEN
        in_data = .TRUE.
        CYCLE
      END IF
      
      IF(.NOT. in_data) CYCLE
      IF(LEN_TRIM(line) == 0) CYCLE
      IF(line(1:1) == '!') CYCLE
      
      ! First non-comment line after header is object name
      IF(LEN_TRIM(name) == 0) THEN
        name = TRIM(ADJUSTL(line))
        CYCLE
      END IF
      
      ! Read Keplerian elements (KEP format line)
      IF(INDEX(line, 'KEP') > 0) THEN
        READ(line,*,IOSTAT=ios) a, e, i, Omega, omega, M, epoch
        IF(ios == 0) EXIT
      END IF
    END DO
    
    CLOSE(unit)
    
    ! Convert angles from degrees to radians
    i = i * DEG2RAD
    Omega = Omega * DEG2RAD
    omega = omega * DEG2RAD
    M = M * DEG2RAD
    
  END SUBROUTINE read_eq1_simple
  
  ! Convert Keplerian elements to Cartesian state
  SUBROUTINE keplerian_to_cartesian(a, e, i, Omega, omega, M, mu, state)
    DOUBLE PRECISION, INTENT(IN) :: a, e, i, Omega, omega, M, mu
    DOUBLE PRECISION, DIMENSION(6), INTENT(OUT) :: state
    
    DOUBLE PRECISION :: E, nu, r, p
    DOUBLE PRECISION :: cos_nu, sin_nu, cos_i, sin_i, cos_Om, sin_Om, cos_om, sin_om
    DOUBLE PRECISION :: x_orb, y_orb, vx_orb, vy_orb
    DOUBLE PRECISION :: n, sqrt_mu_p
    
    ! Solve Kepler's equation for eccentric anomaly
    E = solve_kepler(M, e)
    
    ! True anomaly
    nu = 2.d0 * ATAN2(SQRT(1.d0 + e) * SIN(E/2.d0), SQRT(1.d0 - e) * COS(E/2.d0))
    
    ! Radius
    r = a * (1.d0 - e * COS(E))
    
    ! Position in orbital plane
    cos_nu = COS(nu)
    sin_nu = SIN(nu)
    x_orb = r * cos_nu
    y_orb = r * sin_nu
    
    ! Velocity in orbital plane
    p = a * (1.d0 - e**2)
    sqrt_mu_p = SQRT(mu / p)
    vx_orb = -sqrt_mu_p * sin_nu
    vy_orb = sqrt_mu_p * (e + cos_nu)
    
    ! Rotation to inertial frame
    cos_i = COS(i)
    sin_i = SIN(i)
    cos_Om = COS(Omega)
    sin_Om = SIN(Omega)
    cos_om = COS(omega)
    sin_om = SIN(omega)
    
    ! Position
    state(1) = (cos_Om * cos_om - sin_Om * sin_om * cos_i) * x_orb + &
               (-cos_Om * sin_om - sin_Om * cos_om * cos_i) * y_orb
    state(2) = (sin_Om * cos_om + cos_Om * sin_om * cos_i) * x_orb + &
               (-sin_Om * sin_om + cos_Om * cos_om * cos_i) * y_orb
    state(3) = sin_i * sin_om * x_orb + sin_i * cos_om * y_orb
    
    ! Velocity
    state(4) = (cos_Om * cos_om - sin_Om * sin_om * cos_i) * vx_orb + &
               (-cos_Om * sin_om - sin_Om * cos_om * cos_i) * vy_orb
    state(5) = (sin_Om * cos_om + cos_Om * sin_om * cos_i) * vx_orb + &
               (-sin_Om * sin_om + cos_Om * cos_om * cos_i) * vy_orb
    state(6) = sin_i * sin_om * vx_orb + sin_i * cos_om * vy_orb
    
  END SUBROUTINE keplerian_to_cartesian
  
  ! Solve Kepler's equation: M = E - e*sin(E)
  FUNCTION solve_kepler(M, e) RESULT(E)
    DOUBLE PRECISION, INTENT(IN) :: M, e
    DOUBLE PRECISION :: E
    DOUBLE PRECISION :: E_old, f, fp
    INTEGER :: iter
    INTEGER, PARAMETER :: max_iter = 100
    DOUBLE PRECISION, PARAMETER :: tol = 1.d-12
    
    ! Initial guess
    E = M + e * SIN(M)
    
    ! Newton-Raphson iteration
    DO iter = 1, max_iter
      E_old = E
      f = E - e * SIN(E) - M
      fp = 1.d0 - e * COS(E)
      E = E - f / fp
      
      IF(ABS(E - E_old) < tol) EXIT
    END DO
    
  END FUNCTION solve_kepler
  
  ! Simple 2-body propagation (Keplerian orbit)
  SUBROUTINE propagate_twobody(state, t0, t, mu)
    DOUBLE PRECISION, DIMENSION(6), INTENT(INOUT) :: state
    DOUBLE PRECISION, INTENT(IN) :: t0, t, mu
    
    DOUBLE PRECISION :: dt
    DOUBLE PRECISION :: a, e, i, Omega, omega, M0, M
    
    dt = t - t0
    
    ! Convert state to Keplerian
    CALL cartesian_to_keplerian(state, mu, a, e, i, Omega, omega, M0)
    
    ! Propagate mean anomaly
    M = M0 + SQRT(mu / a**3) * dt
    
    ! Normalize M to [0, 2π]
    M = MOD(M, 2.d0 * PI)
    IF(M < 0.d0) M = M + 2.d0 * PI
    
    ! Convert back to Cartesian
    CALL keplerian_to_cartesian(a, e, i, Omega, omega, M, mu, state)
    
  END SUBROUTINE propagate_twobody
  
  ! Convert Cartesian state to Keplerian elements
  SUBROUTINE cartesian_to_keplerian(state, mu, a, e, i, Omega, omega, M)
    DOUBLE PRECISION, DIMENSION(6), INTENT(IN) :: state
    DOUBLE PRECISION, INTENT(IN) :: mu
    DOUBLE PRECISION, INTENT(OUT) :: a, e, i, Omega, omega, M
    
    DOUBLE PRECISION, DIMENSION(3) :: r_vec, v_vec, h_vec, n_vec, e_vec
    DOUBLE PRECISION :: r, v2, h, n, energy, E, nu
    
    r_vec = state(1:3)
    v_vec = state(4:6)
    
    r = SQRT(DOT_PRODUCT(r_vec, r_vec))
    v2 = DOT_PRODUCT(v_vec, v_vec)
    
    ! Angular momentum
    h_vec(1) = r_vec(2) * v_vec(3) - r_vec(3) * v_vec(2)
    h_vec(2) = r_vec(3) * v_vec(1) - r_vec(1) * v_vec(3)
    h_vec(3) = r_vec(1) * v_vec(2) - r_vec(2) * v_vec(1)
    h = SQRT(DOT_PRODUCT(h_vec, h_vec))
    
    ! Node vector
    n_vec(1) = -h_vec(2)
    n_vec(2) = h_vec(1)
    n_vec(3) = 0.d0
    n = SQRT(DOT_PRODUCT(n_vec, n_vec))
    
    ! Eccentricity vector
    e_vec = ((v2 - mu/r) * r_vec - DOT_PRODUCT(r_vec, v_vec) * v_vec) / mu
    e = SQRT(DOT_PRODUCT(e_vec, e_vec))
    
    ! Semi-major axis
    energy = v2/2.d0 - mu/r
    a = -mu / (2.d0 * energy)
    
    ! Inclination
    i = ACOS(h_vec(3) / h)
    
    ! Longitude of ascending node
    IF(n > 1.d-10) THEN
      Omega = ACOS(n_vec(1) / n)
      IF(n_vec(2) < 0.d0) Omega = 2.d0 * PI - Omega
    ELSE
      Omega = 0.d0
    END IF
    
    ! Argument of perihelion
    IF(n > 1.d-10 .AND. e > 1.d-10) THEN
      omega = ACOS(DOT_PRODUCT(n_vec, e_vec) / (n * e))
      IF(e_vec(3) < 0.d0) omega = 2.d0 * PI - omega
    ELSE
      omega = 0.d0
    END IF
    
    ! True anomaly
    IF(e > 1.d-10) THEN
      nu = ACOS(DOT_PRODUCT(e_vec, r_vec) / (e * r))
      IF(DOT_PRODUCT(r_vec, v_vec) < 0.d0) nu = 2.d0 * PI - nu
    ELSE
      nu = 0.d0
    END IF
    
    ! Eccentric anomaly and mean anomaly
    E = 2.d0 * ATAN(SQRT((1.d0 - e)/(1.d0 + e)) * TAN(nu/2.d0))
    M = E - e * SIN(E)
    
    ! Normalize M
    M = MOD(M, 2.d0 * PI)
    IF(M < 0.d0) M = M + 2.d0 * PI
    
  END SUBROUTINE cartesian_to_keplerian

END PROGRAM orbfit_propagator_simple
