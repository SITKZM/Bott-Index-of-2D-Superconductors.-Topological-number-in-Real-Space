! gfortran -o Bott_index.out Bott_index.f90 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
program Bott_index
    use,intrinsic :: iso_fortran_env
    implicit none
    !parameter
    integer, parameter :: N = 3571, allhop = 7186
    double precision, parameter :: U = 3.1, mu = -1.55, hsq = 0.35, V = 2.55, pi = 4*atan(1.)
    character(7), parameter :: hop_file = "hop.txt"
    character(19), parameter :: Delta_file = "pair_potential.txt", PN_file = "particle_number.txt"
    character(24), parameter :: pos_file = "rescaled_coordinates.txt", E_file = "eigval.txt"
    character(10), parameter :: result_file = "result.txt"
    !variable
    integer :: unit_write_result
    double precision :: h_z = sqrt(hsq)
    double precision :: PN(2*N)
    complex(kind(0d0)) :: Delta(N), Hamiltonian(4 * N, 4 * N)
    double precision :: Bott
    ! for ZHEEVD of lapack
    integer :: INFO, IWORK(3 + 5 * 4 * N)
    double precision :: W(4 * N), RWORK(1 + 5 * 4 * N + 2 * 16 * N**2)
    complex(kind(0d0)) :: WORK(2 * 4 * N + 16 * N**2)
    ! to clock time
    integer(int64) :: time_begin_c,time_end_c, CountPerSec, CountMax

    call system_clock(time_begin_c, CountPerSec, CountMax)

    call read_quantities()

    call make_Hamiltonian()
    call ZHEEVD('V', 'U', 4 * N, Hamiltonian, 4 * N, W, WORK, 2 * 4 * N + 16 * N**2,&
    RWORK, 1 + 20 * N + 2* 16 * N**2, IWORK, 3 + 5 * 4 * N, INFO)
    call get_Bott_index()

    call system_clock(time_end_c)

    call write_files()
contains
    subroutine make_Hamiltonian()
        integer :: i, j, k, l, m
        double precision :: x, y, r
        complex(kind(0d0)) :: H_ij

        Hamiltonian = 0

        ! diagonal elements
        do k = 1, N
            ! spin up
            i = 2 * k - 1
            j = i
            H_ij = -mu - U * PN(2 * k) + h_z

            Hamiltonian(i, j) = H_ij
            Hamiltonian(i + 2 * N, j + 2 * N) = -H_ij

            ! spin down
            i = 2 * k
            j = i
            H_ij = -mu - U * PN(2 * k - 1) - h_z !

            Hamiltonian(i, j) = H_ij
            Hamiltonian(i + 2 * N, j + 2 * N) = -H_ij
        end do

        ! hopping elements
        open(newunit = unit_write_result, file = hop_file)
        do k = 1, allhop
            read(unit_write_result, *) l, m
            l = l + 1
            m = m + 1

            ! spin up
            i = 2 * l - 1
            j = 2 * m - 1
            H_ij = -1.

            Hamiltonian(i, j) = H_ij
            Hamiltonian(j, i) = H_ij
            Hamiltonian(i + 2 * N, j + 2 * N) = -H_ij
            Hamiltonian(j + 2 * N, i + 2 * N) = -H_ij

            ! spin down
            i = 2 * l
            j = 2 * m
            H_ij = -1.

            Hamiltonian(i, j) = H_ij
            Hamiltonian(j, i) = H_ij
            Hamiltonian(i + 2 * N, j + 2 * N) = -H_ij
            Hamiltonian(j + 2 * N, i + 2 * N) = -H_ij
        end do
        close (unit_write_result)

        !spin-orbit coupling
        open (newunit = unit_write_result, file=hop_file)
        do k = 1, allhop
            read (unit_write_result, *) l, m, x, y
            l = l + 1
            m = m + 1

            r = sqrt(x**2 + y**2)
            x = x/r
            y = y/r

            !iup, jdown
            i = 2*l - 1
            j = 2*m
            H_ij = V*cmplx(-x, y, kind(0d0))

            Hamiltonian(i, j) = H_ij
            Hamiltonian(j, i) = conjg(H_ij)
            Hamiltonian(i + 2 * N, j + 2 * N) = -conjg(H_ij)
            Hamiltonian(j + 2 * N, i + 2 * N) = -H_ij

            !idown, jup
            i = 2*l
            j = 2*m - 1
            H_ij = V*cmplx(x, y, kind(0d0))

            Hamiltonian(i, j) = H_ij
            Hamiltonian(j, i) = conjg(H_ij)
            Hamiltonian(i + 2 * N, j + 2 * N) = -conjg(H_ij)
            Hamiltonian(j + 2 * N, i + 2 * N) = -H_ij
        end do
        close(unit_write_result)

        !pair potential
        do k = 1, N
            !(H_ij)c_{idown}^dag c_{iup}^dag
            i = 2*k
            j = 2*k - 1 + 2*N
            H_ij = -Delta(k)

            Hamiltonian(i, j) = H_ij
            Hamiltonian(j, i) = conjg(H_ij)

            !(H_ij)c_{iup}^dag c_{idown}^dag
            i = 2*k - 1
            j = 2*k + 2*N
            H_ij = Delta(k)

            Hamiltonian(i, j) = H_ij
            Hamiltonian(j, i) = conjg(H_ij)
        end do
    end subroutine make_Hamiltonian

    subroutine get_Bott_index()
        integer :: i, VL, VR
        double precision :: pos_x, pos_y, RWORK_2(2 * 2 * N)
        complex(kind(0d0)) :: X(4 * N, 4 * N), Y(4 * N, 4 * N)
        complex(kind(0d0)) :: P1(2 * N, 4 * N), P2(4 * N, 2 * N)
        complex(kind(0d0)) :: U_X(2 * N, 2 * N), U_Y(2 * N, 2 * N), U_prod(2 * N, 2 * N)
        complex(kind(0d0)) :: W_2(2 * N)

        X = 0
        Y = 0
        P1 = 0
        P2 = 0

        !position operators
        open(newunit=unit_write_result, file=pos_file)
        do i = 1, N
            read(unit_write_result, *) pos_x, pos_y
            X(2 * i - 1, 2 * i - 1) = exp(cmplx(0, 2 * pi * pos_x, kind(0d0)))
            X(2 * i, 2 * i) = exp(cmplx(0, 2 * pi * pos_x, kind(0d0)))
            X(2 * i - 1 + 2 * N, 2 * i - 1 + 2 * N) = exp(cmplx(0, 2 * pi * pos_x, kind(0d0)))
            X(2 * i + 2 * N, 2 * i + 2 * N) = exp(cmplx(0, 2 * pi * pos_x, kind(0d0)))

            Y(2 * i - 1, 2 * i - 1) = exp(cmplx(0, 2 * pi * pos_y, kind(0d0)))
            Y(2 * i, 2 * i) = exp(cmplx(0, 2 * pi * pos_y, kind(0d0)))
            Y(2 * i - 1 + 2 * N, 2 * i - 1 + 2 * N) = exp(cmplx(0, 2 * pi * pos_y, kind(0d0)))
            Y(2 * i + 2 * N, 2 * i + 2 * N) = exp(cmplx(0, 2 * pi * pos_y, kind(0d0)))
        end do
        close(unit_write_result)

        do i = 1, 2 * N
            P1(i, i) = 1
            P2(i, i) = 1
        end do

        U_X = matmul(P1, matmul(transpose(conjg(Hamiltonian)), matmul(X, matmul(Hamiltonian, P2))))
        U_Y = matmul(P1, matmul(transpose(conjg(Hamiltonian)), matmul(Y, matmul(Hamiltonian, P2))))
        U_prod = matmul(U_Y, matmul(U_X, matmul(transpose(conjg(U_Y)), transpose(conjg(U_X)))))

        !Bott index
        Bott = 0
        call ZGEEV('N', 'N', 2 * N, U_prod, 2 * N, W_2, VL, 1, VR, 1, WORK, 2*2*N, RWORK_2, INFO)
        do i = 1, 2*N
            Bott = Bott + aimag(log(W_2(i)/abs(W_2(i))))
        end do
        Bott = Bott / (2 * pi)
    end subroutine get_Bott_index

    subroutine write_files()
        integer :: i

        ! eigenenergies
        open(newunit=unit_write_result, file=E_file)
            do i = 1, 4 * N
                write(unit_write_result, *) W(i)
            end do
        close(unit_write_result)

        ! hsq max_Delta mean_Delta min_Delta min_E Bott
        open(newunit=unit_write_result, file=result_file)
            write(unit_write_result, *) "Run time=", real(time_end_c - time_begin_c, kind(0d0)) / CountPerSec, "sec"
            write(unit_write_result, *) "hsq, max_Delta, mean_Delta, min_Delta, min_E, Bott:"
            write(unit_write_result, *) hsq
            write(unit_write_result, *) maxval(abs(Delta))
            write(unit_write_result, *) sum(abs(Delta)) / N
            write(unit_write_result, *) minval(abs(Delta))
            write(unit_write_result, *) minval(abs(W))
            write(unit_write_result, *) Bott
        close(unit_write_result)
    end subroutine write_files

    subroutine read_quantities()
        integer :: i
        double precision :: s, t

        open(newunit = unit_write_result, file = Delta_file)
        do i = 1, N
            read(unit_write_result, *) s, t
            Delta(i) = cmplx(s * cos(t), s * sin(t), kind(0d0))
        end do
        close(unit_write_result)

        open(newunit = unit_write_result, file = PN_file)
        do i = 1, N
            read(unit_write_result, *) s, t
            PN(2 * i - 1) = s
            PN(2 * i) = t
        end do
        close(unit_write_result)
    end subroutine read_quantities
end program Bott_index
