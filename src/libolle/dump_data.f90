#include "precompilerdefinitions"
module dump_data
!! Routines to dump data to files, and produce a .gnuplot file that can plot the data.
use konstanter, only: flyt,lo_gnuplotterminal
use gottochblandat, only: open_file,tochar,qsort
implicit none

private
public :: lo_dump_gnuplot_2d_real
public :: lo_dump_palette_to_gnuplot

contains

!> Dumps a clean data file with an associated gnuplot command file. If add_x is set to 'yes' an extra column with x-values is added. @todo Perhaps clean this a little, maybe more useful if the input is x(:),y(:,:). Also, it could be nice to scale the output to some unit without touching the original data.
subroutine lo_dump_gnuplot_2d_real(m,filename,add_x,xlabel,ylabel,title,legend,palette,yformat,xformat,linestyle,xrange,yrange,keyplacement)
    !> Data to be plotted
    real(flyt), dimension(:,:), intent(in) :: m
    !> Filename. The gnuplot file will have the same name with suffix .gnuplot
    character(len=*), intent(in) :: filename
    !> should an x-axis be added?
    logical, intent(in), optional :: add_x
    !> Label on x-axis
    character(len=*), intent(in), optional :: xlabel
    !> Label on y-axis
    character(len=*), intent(in), optional :: ylabel
    !> Title of plot
    character(len=*), intent(in), optional :: title
    !> Line colors. Choose bluegold or redgreen for smooth colors, and default for high contrast lines.
    character(len=*), intent(in), optional :: palette
    !> gnuplot syntax format for x-axis
    character(len=*), intent(in), optional :: xformat
    !> gnuplot syntax format for y-axis
    character(len=*), intent(in), optional :: yformat
    !> linestyle in gnuplot syntax
    character(len=*), intent(in), optional :: linestyle
    !> where the legend should be, gnuplot syntax
    character(len=*), intent(in), optional :: keyplacement
    !> min and max for the x-axis, gnuplot syntax
    real(flyt), dimension(2), intent(in), optional :: xrange
    !> min and max for y-axis, gnuplot syntax
    real(flyt), dimension(2), intent(in), optional :: yrange
    !> legend, should have one entry per column in the data
    character(len=*), dimension(:), intent(in), optional :: legend

    integer :: i,j,u,ncol
    character(len=1000) :: opf,command_filename,ls

    ncol=size(m,1)
    if ( present(add_x) ) then
        if ( add_x ) ncol=ncol+1
    endif

    command_filename=trim(filename)//'.gnuplot'
    u=open_file('out',command_filename)
        ! gnuplot defaults to 90 dpi, 355 pixels is then 10cm, a nice size.
        ! enhanced makes some formatting syntax possible, e.g. 3^2 = 3²
        write(u,*) 'set terminal '//trim(lo_gnuplotterminal)//' size 500,350 enhanced font "CMU Serif,10"'
        write(u,*) 'set border lw 0.5'
        if ( present(keyplacement) ) then
            write(u,*) trim(keyplacement)
        else
            write(u,*) 'set key inside left Left reverse'
        endif
        write(u,*) '#set rmargin 4'
        write(u,*) '#set lmargin 8'
        write(u,*) '#set bmargin 3'
        write(u,*) '#set xlabel offset 0,0.35'
        write(u,*) '#set ylabel offset 0.5,0'
        write(u,*) 'set ytics scale 0.5'
        write(u,*) 'set xtics scale 0.5'
        write(u,*) 'set mytics 10'
        write(u,*) 'set mxtics 10'

        ! Maybe format the axes
        if ( present(yformat) ) then
            if ( trim(yformat) .eq. 'none' ) then
                write(u,*) 'unset ytics'
            else
                write(u,*) 'set format y '//trim(yformat)
            endif
        endif
        if ( present(xformat) ) then
            if ( trim(xformat) .eq. 'none' ) then
                write(u,*) 'unset xtics'
            else
                write(u,*) 'set format x '//trim(yformat)
            endif
        endif

        ! maybe ranges
        if ( present(xrange) ) write(u,*) 'set xrange [',real(xrange(1)),':',real(xrange(2)),']'
        if ( present(yrange) ) write(u,*) 'set yrange [',real(yrange(1)),':',real(yrange(2)),']'

        ! Maybe nonstandard linestyle
        if ( present(linestyle) ) then
            ls=trim(linestyle)
        else
            ls='line'
        endif

        ! Set my nice palettes, maybe
        if ( present(palette) ) then
            select case(trim(palette))
                case('bluegold')
                    write(u,*) 'set style line 1 lc rgb "#917318" '
                    write(u,*) 'set style line 2 lc rgb "#B4901A" '
                    write(u,*) 'set style line 3 lc rgb "#D9C86A" '
                    write(u,*) 'set style line 4 lc rgb "#E9DEA1" '
                    write(u,*) 'set style line 5 lc rgb "#B3E4CD" '
                    write(u,*) 'set style line 6 lc rgb "#7FD3C0" '
                    write(u,*) 'set style line 7 lc rgb "#64CAB4" '
                    write(u,*) 'set style line 8 lc rgb "#47B8BA" '
                    write(u,*) 'set style line 9 lc rgb "#3788AB" '
                    write(u,*) 'set style line 10 lc rgb "#286388"'
                    write(u,*) 'set style increment user'
                case('redgreen')
                    write(u,*) 'set style line 1 lc rgb "#4A141E" '
                    write(u,*) 'set style line 2 lc rgb "#88172B" '
                    write(u,*) 'set style line 3 lc rgb "#AA253A" '
                    write(u,*) 'set style line 4 lc rgb "#B4616F" '
                    write(u,*) 'set style line 5 lc rgb "#DFB1B0" '
                    write(u,*) 'set style line 6 lc rgb "#C3C055" '
                    write(u,*) 'set style line 7 lc rgb "#90AB3F" '
                    write(u,*) 'set style line 8 lc rgb "#317E41" '
                    write(u,*) 'set style line 9 lc rgb "#1F6341" '
                    write(u,*) 'set style line 10 lc rgb "#143E29" '
                    write(u,*) 'set style increment user'
            end select
        else ! use my default 9 high-contrast colors
            write(u,*) 'set style line 1 lc rgb "#e41a1c" lw 1'
            write(u,*) 'set style line 2 lc rgb "#377eb8" lw 1'
            write(u,*) 'set style line 3 lc rgb "#4daf4a" lw 1'
            write(u,*) 'set style line 4 lc rgb "#984ea3" lw 1'
            write(u,*) 'set style line 5 lc rgb "#ff7f00" lw 1'
            write(u,*) 'set style line 6 lc rgb "#ffff33" lw 1'
            write(u,*) 'set style line 7 lc rgb "#a65628" lw 1'
            write(u,*) 'set style line 8 lc rgb "#f781bf" lw 1'
            write(u,*) 'set style line 9 lc rgb "#999999" lw 1'
            write(u,*) 'set style increment user'
        endif
        !
        if ( present(title) ) write(u,*) 'set title "',trim(title),'" '
        if ( present(xlabel) ) write(u,*) 'set xlabel "',trim(xlabel),'" '
        if ( present(ylabel) ) write(u,*) 'set ylabel "',trim(ylabel),'" '
        write(u,'(A)',advance='no') 'plot'
        do j=2,ncol
            write(u,'(A)',advance='no') ' "'//trim(filename)//'" u 1:'//tochar(j)//' w '//trim(ls)
            if( present(legend) ) then
                write(u,'(A)',advance='no') ' t "'//trim(adjustl(legend(j-1)))//'"'
            endif
            if( j .lt. ncol ) then
                write(u,'(A)',advance='no') ","
            endif
        enddo
        write(u,*) ' '
        !
    close(u)

    ! write the actual data
    opf="(1X,"//tochar(ncol)//"(1X,E20.10))"
    u=open_file('out',trim(filename))
        do i=1,size(m,2)
            if ( present(add_x) ) then
                if ( add_x ) then
                    write(u,opf) i*1.0_flyt,m(:,i)
                endif
            else
                write(u,opf) m(:,i)
            endif
        enddo
    close(u)
end subroutine lo_dump_gnuplot_2d_real

!> Writes nice color scale in gnuplot format Adaptation of http://www.mathworks.com/matlabcentral/fileexchange/25690-haxby-color-map. Colormap is based on the colors used by W. F. Haxby's Gravity field of World's oceans, 1985, developed for geoid and gravity maps. The version used here is formed from a linear interpolation of the GMT color table used by MB-System by David W. Caress and Dale N. Chayes (http://www.ldeo.columbia.edu/res/pi/MB-System). Writes the palette to gnuplot file, passed as the unit the gnuplot file is opened with.
subroutine lo_dump_palette_to_gnuplot(u,logscale)
    !> Output unit for gnuplot file
    integer, intent(in) :: u
    !> Should the color scale be logarithmic?
    logical, intent(in), optional :: logscale
    !
    character(len=20), dimension(100) :: clr
    real(flyt), dimension(100) :: dum
    real(flyt) :: f0,f1,f2
    logical :: ls
    integer :: i

    if ( present(logscale) ) then
        ls=logscale
    else
        ls=.false.
    endif

    clr(  1)='"#440154"'
    clr(  2)='"#450458"'
    clr(  3)='"#46085C"'
    clr(  4)='"#460C5F"'
    clr(  5)='"#471063"'
    clr(  6)='"#471466"'
    clr(  7)='"#48176A"'
    clr(  8)='"#481B6D"'
    clr(  9)='"#481E70"'
    clr( 10)='"#482173"'
    clr( 11)='"#482575"'
    clr( 12)='"#472878"'
    clr( 13)='"#472B7A"'
    clr( 14)='"#472E7C"'
    clr( 15)='"#46317F"'
    clr( 16)='"#453580"'
    clr( 17)='"#443882"'
    clr( 18)='"#433B84"'
    clr( 19)='"#423E85"'
    clr( 20)='"#414186"'
    clr( 21)='"#404487"'
    clr( 22)='"#3F4788"'
    clr( 23)='"#3E4A89"'
    clr( 24)='"#3D4D8A"'
    clr( 25)='"#3B508B"'
    clr( 26)='"#3A528B"'
    clr( 27)='"#39558C"'
    clr( 28)='"#38588C"'
    clr( 29)='"#365B8D"'
    clr( 30)='"#355D8D"'
    clr( 31)='"#34608D"'
    clr( 32)='"#33638D"'
    clr( 33)='"#31658E"'
    clr( 34)='"#30688E"'
    clr( 35)='"#2F6A8E"'
    clr( 36)='"#2E6D8E"'
    clr( 37)='"#2D6F8E"'
    clr( 38)='"#2C728E"'
    clr( 39)='"#2B748E"'
    clr( 40)='"#2A778E"'
    clr( 41)='"#29798E"'
    clr( 42)='"#287C8E"'
    clr( 43)='"#277E8E"'
    clr( 44)='"#26808E"'
    clr( 45)='"#25838E"'
    clr( 46)='"#24858E"'
    clr( 47)='"#23888E"'
    clr( 48)='"#228A8D"'
    clr( 49)='"#228D8D"'
    clr( 50)='"#218F8D"'
    clr( 51)='"#20918C"'
    clr( 52)='"#1F948C"'
    clr( 53)='"#1F968B"'
    clr( 54)='"#1E998B"'
    clr( 55)='"#1E9B8A"'
    clr( 56)='"#1E9E89"'
    clr( 57)='"#1EA088"'
    clr( 58)='"#1FA287"'
    clr( 59)='"#20A586"'
    clr( 60)='"#21A785"'
    clr( 61)='"#23AA83"'
    clr( 62)='"#25AC82"'
    clr( 63)='"#28AE80"'
    clr( 64)='"#2AB17E"'
    clr( 65)='"#2EB37D"'
    clr( 66)='"#31B57B"'
    clr( 67)='"#35B779"'
    clr( 68)='"#39BA76"'
    clr( 69)='"#3DBC74"'
    clr( 70)='"#42BE72"'
    clr( 71)='"#46C06F"'
    clr( 72)='"#4BC26C"'
    clr( 73)='"#50C469"'
    clr( 74)='"#56C666"'
    clr( 75)='"#5BC863"'
    clr( 76)='"#61CA60"'
    clr( 77)='"#66CC5D"'
    clr( 78)='"#6CCE59"'
    clr( 79)='"#72D055"'
    clr( 80)='"#78D152"'
    clr( 81)='"#7FD34E"'
    clr( 82)='"#85D54A"'
    clr( 83)='"#8CD646"'
    clr( 84)='"#92D841"'
    clr( 85)='"#99D93D"'
    clr( 86)='"#A0DA39"'
    clr( 87)='"#A7DB34"'
    clr( 88)='"#ADDD30"'
    clr( 89)='"#B4DE2B"'
    clr( 90)='"#BBDF27"'
    clr( 91)='"#C2E023"'
    clr( 92)='"#C9E11F"'
    clr( 93)='"#D0E21C"'
    clr( 94)='"#D7E319"'
    clr( 95)='"#DEE318"'
    clr( 96)='"#E4E418"'
    clr( 97)='"#EBE51A"'
    clr( 98)='"#F1E61C"'
    clr( 99)='"#F8E720"'
    clr(100)='"#FEE724"'

    if ( ls ) then
        f1=1E-2_flyt
        do i=1,100
            dum(i)=log( (i-1.0_flyt)/100.0_flyt + f1 )
        enddo
        dum=dum-minval(dum)
        dum=dum/maxval(dum)
        dum=1.0_flyt-dum
        call qsort(dum)

        f2=log(f1)
        write(u,'(A)',advance='no') 'set palette defined ('
        do i=1,99
            f0=dum(i)
            write(u,'(F18.12," ",A,",")',advance='no') f0,trim(clr(i))
        enddo
        write(u,'(F18.12," ",A,")")') 1.0_flyt,trim(clr(i))
    else
        write(u,'(A)',advance='no') 'set palette defined ('
        do i=1,99
            write(u,'(F18.12," ",A,",")',advance='no') (i)*1.0_flyt,trim(clr(i))
        enddo
        write(u,'(F18.12," ",A,")")') (100)*1.0_flyt,trim(clr(i))
    endif
end subroutine lo_dump_palette_to_gnuplot

end module
