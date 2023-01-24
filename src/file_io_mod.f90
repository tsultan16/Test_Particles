MODULE file_io_mod

USE constants_mod
USE grid_arrays_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE file_output(t_dump)

    INTEGER, INTENT(IN) :: t_dump
    INTEGER :: i
    CHARACTER (LEN=200) :: filename
    CHARACTER(LEN=6) :: uniti


    IF(t_dump<10) THEN
        WRITE(uniti,'(I1.1)') t_dump
    ELSE IF(t_dump>=10 .and. t_dump<100) THEN
        WRITE(uniti,'(I2.2)') t_dump
    ELSE IF(t_dump>=100 .and. t_dump<1000) THEN
        WRITE (uniti,'(I3.3)') t_dump
    ELSE IF(t_dump>=1000 .and. t_dump<10000) THEN
        WRITE (uniti,'(I4.3)') t_dump
    ELSE IF(t_dump>=10000 .and. t_dump<100000) THEN
        WRITE (uniti,'(I5.3)') t_dump  
    END IF   

        
    filename = TRIM('Output/particles_dump=')//TRIM(uniti)//TRIM('.dat')
    OPEN(UNIT = 1, FILE = filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM', STATUS = 'NEW')
    
    ! write the electron data first
    DO i = 1, Ne
        WRITE(1) electron(i)%x,electron(i)%v 
    END DO
        
    ! now write the protons
    DO i = 1, Np
        WRITE(1) proton(i)%x,proton(i)%v 
    END DO
    
    CLOSE(UNIT=1)


END SUBROUTINE file_output






END MODULE file_io_mod