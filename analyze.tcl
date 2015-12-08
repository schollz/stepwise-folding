mol new model1_protein_autopsf.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile eq1.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

set mol 0
set file [open "rmsd.dat" w]         
# use frame 0 for the reference
set reference [atomselect $mol "protein" frame 0]
# the frame being compared
set compare [atomselect $mol "protein"]

set num_steps [molinfo $mol get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
        # get the correct frame
        $compare frame $frame

        # compute the transformation
        set trans_mat [measure fit $compare $reference]
        # do the alignment
        $compare move $trans_mat
        # compute the RMSD
        set rmsd [measure rmsd $compare $reference]
        # print the RMSD
        puts $file "$frame $rmsd"
}

exit