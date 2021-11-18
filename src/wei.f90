module wei

    interface 
        function threej(j1,j2,j3,m1,m2,m3,pprint)
            implicit none
            real j1,j2,j3,m1,m2,m3
            logical pprint
            real threej
        end function threej

        function sixj(a,b,c,d,e,f,pprint)
            implicit none          
            real a,b,c,d,e,f
            logical pprint
            real sixj
        end function sixj

        function ninej(a,b,c,d,e,f,g,h,j,pprint)
            implicit none
            real a,b,c,d,e,f,g,h,j
            logical pprint
            real ninej
        end function ninej
    end interface
end module wei        
