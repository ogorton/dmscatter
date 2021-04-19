program darkmatter

    use main

    implicit none

    integer :: computeoption

    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'Welcome to our FORTAN 90 implementation of the WIMP-nucleon scattering code.'
    print*,'Based on the Mathematica script described in arXiv:1308.6288 (2003).'
    print*,'  VERSION 1.4 UPDATED: 2021.01.15 @ SDSU'
    print*,'  Dev. contact: cjohnson@sdsu.edu'
    print*,'                 ogorton@sdsu.edu'
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,''
    print*,'Select an option:'
    print*,'[1] Event rate per unit recoil energy (spectra)'
    print*,'[2] Scattering probability'
    print*,'[3] Differential cross section per recoil energy'
    print*,'[4] (*) Total cross section'
    print*,'[5] Total scattering events per detector (does not produce spectra data)'
    print*," "
    print*,"(*) Not available in this release."
    read*,computeoption

    call computefunction(computeoption)

end program  

