
proc BuildGraphene {lx ly bl boxz type nlayer isbond} {

    graphene -lx $lx -ly $ly -type $type -nlayers $nlayer -cc $bl -b $isbond -a $isbond -d $isbond -i $isbond -t C-C

    set all [atomselect top all]

    set xmin [lindex [measure minmax $all] 0 0]
    set ymin [lindex [measure minmax $all] 0 1]
    $all moveby "[expr -$xmin] [expr -$ymin] 0.0"

    set xmin [lindex [measure minmax $all] 0 0]
    set ymin [lindex [measure minmax $all] 0 1]
    set xmax [lindex [measure minmax $all] 1 0]
    set ymax [lindex [measure minmax $all] 1 1]

    if { $nlayer == 1 } {

        if { $type == "armchair" } {
            set boxx [expr $xmax-$xmin+$bl*5*sqrt(3.0)]
            set boxy [expr $ymax-$ymin+$bl*10]
        } elseif { $type == "zigzag" } {
            set boxx [expr $xmax-$xmin+$bl*10]
            set boxy [expr $ymax-$ymin+$bl*5*sqrt(3.0)]
        }

    } else {

        for { set a 0 }  {$a < $nlayer} {incr a} {

            set nlaytop [expr $a * 3.35 + 0.1]
            set nlaybot [expr $a * 3.35 - 0.1]
            set newtype [expr $a + 1]
            [atomselect top "z>$nlaybot and z<$nlaytop"] set name $newtype
            [atomselect top "z>$nlaybot and z<$nlaytop"] set type $newtype
            [atomselect top "z>$nlaybot and z<$nlaytop"] set element $newtype

        }

        if { $type == "armchair" } {
            set boxx [expr $xmax-$xmin+$bl*5*sqrt(3.0)]
            set boxy [expr $ymax-$ymin]
        } elseif { $type == "zigzag" } {
            set boxx [expr $xmax-$xmin+$bl*5]
            set boxy [expr $ymax-$ymin]
        }

    }

    pbc set "$boxx $boxy $boxz"

    return [[atomselect top all] molid]

}


# lenth units: nm
set lx 5
set ly 5
set bl 0.142
set type "armchair"
set nlayer 1
set isbond 0

# units: A, default layer distance: 3.35, which cannot be modified.
set boxz [expr $nlayer * 3.35]

# Single layer graphene. You should define the boxz firstly.
# set boxz 10.0
set multigra [BuildGraphene $lx $ly $bl $boxz $type $nlayer $isbond]

topo writelammpsdata read.data full


mol modmaterial 0 $multigra Diffuse
mol material Diffuse
mol modrep 0 $multigra
mol modstyle 0 $multigra VDW 0.200000 50.000000
mol addrep $multigra
mol modstyle 1 $multigra DynamicBonds 1.600000 0.200000 12.000000

display projection perspective
display projection orthographic
display depthcue off

color Display Background white

#topo guessangles
#topo guessdihedrals
#topo guessimpropers
