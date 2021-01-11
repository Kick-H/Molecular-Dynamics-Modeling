# length units: nm
set lx 25
set ly 25
set bl 0.142
set boxz 10.0

set type "armchair"
set type "zigzag"

graphene -lx $lx -ly $ly -type $type -cc $bl -b 0 -a 0 -d 0 -i 0 -t C-C

set all [atomselect top all]
set xmin [lindex [measure minmax $all] 0 0]
set ymin [lindex [measure minmax $all] 0 1]
$all moveby "[expr -$xmin] [expr -$ymin] 0.0"

set xmin [lindex [measure minmax $all] 0 0]
set ymin [lindex [measure minmax $all] 0 1]
set xmax [lindex [measure minmax $all] 1 0]
set ymax [lindex [measure minmax $all] 1 1]
#set xcen [lindex [measure center $all] 0]

if { $type == "armchair" } {
    set boxx [expr $xmax-$xmin+$bl*5*sqrt(3.0)]
    set boxy [expr $ymax-$ymin+$bl*10]
} elseif { $type == "zigzag" } {
    set boxx [expr $xmax-$xmin+$bl*10]
    set boxy [expr $ymax-$ymin+$bl*5*sqrt(3.0)]
}

pbc set "$boxx $boxy $boxz"

topo writelammpsdata gra-$lx-$ly.data atomic

exit
