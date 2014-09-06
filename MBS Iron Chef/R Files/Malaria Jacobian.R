J = matrix( c(-bitel*IL-biteh*IH, vecdeath, vecdeath, 0, 0, -bitel*Sv, -biteh*Sv, 0,
              bitel*IL + bite*IH, -vecdeath - vecexpostime, 0, 0, 0, bitel*Sv, biteh*Sv, 0,
              0, vecexpostime, vecdeath, 0, 0, 0, 0, 0,
              0, 0, -bites*S, -sdeath-bites+births, birth, birth, birth, 0,
              0, 0, bites*S, bites * Iv, -edeath-etime, 0, 0, 0,
              0, 0, 0, 0, etime, -iltime-ildeath, ihtoil, 0,
              0, 0, 0, 0, 0, iltime, -ihdeath-ihtoil, 0,
              0, 0, 0, 0 , 0, ildeath, ihdeath, 0), nrow = 8, byrow = TRUE)