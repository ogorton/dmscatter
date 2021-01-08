function dmresponsefun(eft, term, tau1, tau2, q, v, jchi, muT)
    use kinds
    use parameters
    implicit none
    interface
        function dmresponse_iso(eft, term, tau1, tau2, q, v, jchi, muT)
            use kinds
            use parameters
            real(doublep), allocatable, intent(in) :: eft(:,:)
            integer :: term
            integer :: tau1, tau2
            real(doublep) :: q, v, jchi, muT
            REAL(doublep) :: dmresponse_iso
        end function
        function dmresponse_pn(eft, term, tau1, tau2, q, v, jchi, muT)
            use kinds
            use parameters
            real(doublep), allocatable, intent(in) :: eft(:,:)
            integer :: term
            integer :: tau1, tau2
            real(doublep) :: q, v, jchi, muT
            REAL(doublep) :: dmresponse_pn
        end function        
    end interface

    real(doublep), allocatable, intent(in) :: eft(:,:)
    integer :: term
    integer :: tau1, tau2
    real(doublep) :: q, v, jchi, muT
    REAL(doublep) :: dmresponsefun

    if (pndens) then
        dmresponsefun = dmresponse_pn(eft, term, tau1, tau2, q, v, jchi, muT)
    else
        dmresponsefun = dmresponse_iso(eft, term, tau1, tau2, q, v, jchi, muT)
    end if
    
end function
