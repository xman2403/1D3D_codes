        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 13 16:37:05 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MYTRIM__genmod
          INTERFACE 
            FUNCTION MYTRIM(STR) RESULT(RES)
              CHARACTER(*) ,TARGET :: STR
              CHARACTER(LEN=:) ,POINTER :: RES
            END FUNCTION MYTRIM
          END INTERFACE 
        END MODULE MYTRIM__genmod
