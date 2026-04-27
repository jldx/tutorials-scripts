# Tcl script for VMD: Extract coordinates of a selection across a trajectory for a specific molecule
# Usage: source this script in VMD's TkConsole, then call extract_coordinates

proc extract_coordinates {molid selection axis outfilename {startframe 0} {endframe -1} {stride 1}} {
    # Check if the molecule ID is valid by attempting to get its name
    if {[catch {molinfo $molid get name}]} {
        puts "Error: Molecule ID $molid does not exist."
        return
    }

    # Check if the selection is valid
    set sel [atomselect $molid "$selection" frame $startframe]
    if {[$sel num] == 0} {
        puts "Error: No atoms selected. Check your selection string."
        $sel delete
        return
    }
    $sel delete

    # Set endframe to the last frame if not specified
    set nframes [molinfo $molid get numframes]
    if {$endframe == -1 || $endframe >= $nframes} {
        set endframe [expr {$nframes - 1}]
    }

    # Open output file
    set out [open $outfilename w]

    # Write header
    if {$axis eq "xyz"} {
        puts $out "# Frame\tAtomID\tResID\tResName\tAtomName\tX\tY\tZ"
    } else {
        puts $out "# Frame\tAtomID\tResID\tResName\tAtomName\t$axis"
    }

    # Loop over frames
    for {set frame $startframe} {$frame <= $endframe} {incr frame $stride} {
        set sel [atomselect $molid "$selection" frame $frame]
        set n [$sel num]

        # Get all atom properties at once for efficiency
        set atomIDs [$sel get index]
        set resIDs [$sel get resid]
        set resNames [$sel get resname]
        set atomNames [$sel get name]

        if {$axis eq "xyz"} {
            set xCoords [$sel get x]
            set yCoords [$sel get y]
            set zCoords [$sel get z]
        } else {
            set coords [$sel get $axis]
        }

        # Loop over selected atoms
        for {set i 0} {$i < $n} {incr i} {
            set atomID [lindex $atomIDs $i]
            set resID [lindex $resIDs $i]
            set resName [lindex $resNames $i]
            set atomName [lindex $atomNames $i]

            if {$axis eq "xyz"} {
                set x [lindex $xCoords $i]
                set y [lindex $yCoords $i]
                set z [lindex $zCoords $i]
                puts $out "$frame\t$atomID\t$resID\t$resName\t$atomName\t$x\t$y\t$z"
            } else {
                set coord [lindex $coords $i]
                puts $out "$frame\t$atomID\t$resID\t$resName\t$atomName\t$coord"
            }
        }
        $sel delete
    }

    close $out
    puts "Coordinates extracted to $outfilename"
}

# Example usage:
# extract_coordinates 0 "resid 244 245 246 and name CA" xyz "my_output.dat"
# This will extract xyz-coordinates of CA atoms in residues 244, 245, and 246 for molecule 0 and save to "my_output.dat".