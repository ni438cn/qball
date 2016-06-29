set nrowmax 48
si54p.sys
set ecut 65.0

set wf_dyn PSDA
set ecutprec 6.0
set xc PBE

symmetry  1  0  0  0  1  0  0  0  1
symmetry  0  1  -1  1  0  -1  0  0  -1  0.3333333333  -0.3333333333  0.0000000
symmetry  -1  0  0  -1  0  1  -1  1  0  0.0000000  -0.3333333333  0.3333333333
symmetry  0  -1  1  0  -1  0  1  -1  0  0.0000000  -0.3333333333  0.0000000
symmetry  -1  0  1  -1  1  0  -1  0  0
symmetry  0  1  0  0  0  1  1  0  0  0.0000000  -0.3333333333  0.3333333333
symmetry  1  0  -1  0  0  -1  0  1  -1  0.0000000  -0.3333333333  0.0000000
symmetry  0  -1  0  1  -1  0  0  -1  1  0.3333333333  -0.3333333333  0.0000000
symmetry  0  0  -1  0  1  -1  1  0  -1
symmetry  -1  1  0  -1  0  0  -1  0  1  0.0000000  -0.3333333333  0.0000000
symmetry  0  0  1  1  0  0  0  1  0  0.3333333333  -0.3333333333  0.0000000
symmetry  1  -1  0  0  -1  1  0  -1  0  0.0000000  -0.3333333333  0.3333333333
symmetry  0  1  0  1  0  0  0  0  1  0.3333333333  -0.3333333333  0.0000000
symmetry  1  0  -1  0  1  -1  0  0  -1
symmetry  0  -1  0  0  -1  1  1  -1  0  0.0000000  -0.3333333333  0.3333333333
symmetry  -1  0  1  -1  0  0  -1  1  0  0.0000000  -0.3333333333  0.0000000
symmetry  1  0  0  0  0  1  0  1  0  0.0000000  -0.3333333333  0.3333333333
symmetry  -1  0  0  -1  1  0  -1  0  1
symmetry  0  -1  1  1  -1  0  0  -1  0  0.3333333333  -0.3333333333  0.0000000
symmetry  0  1  -1  0  0  -1  1  0  -1  0.0000000  -0.3333333333  0.0000000
symmetry  1  -1  0  0  -1  0  0  -1  1  0.0000000  -0.3333333333  0.0000000
symmetry  0  0  1  0  1  0  1  0  0
symmetry  0  0  -1  1  0  -1  0  1  -1  0.3333333333  -0.3333333333  0.0000000
symmetry  -1  1  0  -1  0  1  -1  0  0  0.0000000  -0.3333333333  0.3333333333

set nparallelkpts 4
set nkpoints 29
kpoint   0.0000000   0.0000000   0.0000000   0.0039063  crystal
kpoint   0.0000000   0.0000000   0.1250000   0.0312500  crystal
kpoint   0.0000000   0.0000000   0.2500000   0.0312500  crystal
kpoint   0.0000000   0.0000000   0.3750000   0.0312500  crystal
kpoint   0.0000000   0.0000000  -0.5000000   0.0156250  crystal
kpoint   0.0000000   0.1250000   0.1250000   0.0234375  crystal
kpoint   0.0000000   0.1250000   0.2500000   0.0937500  crystal
kpoint   0.0000000   0.1250000   0.3750000   0.0937500  crystal
kpoint   0.0000000   0.1250000  -0.5000000   0.0937500  crystal
kpoint   0.0000000   0.1250000  -0.3750000   0.0937500  crystal
kpoint   0.0000000   0.1250000  -0.2500000   0.0937500  crystal
kpoint   0.0000000   0.1250000  -0.1250000   0.0468750  crystal
kpoint   0.0000000   0.2500000   0.2500000   0.0234375  crystal
kpoint   0.0000000   0.2500000   0.3750000   0.0937500  crystal
kpoint   0.0000000   0.2500000  -0.5000000   0.0937500  crystal
kpoint   0.0000000   0.2500000  -0.3750000   0.0937500  crystal
kpoint   0.0000000   0.2500000  -0.2500000   0.0468750  crystal
kpoint   0.0000000   0.3750000   0.3750000   0.0234375  crystal
kpoint   0.0000000   0.3750000  -0.5000000   0.0937500  crystal
kpoint   0.0000000   0.3750000  -0.3750000   0.0468750  crystal
kpoint   0.0000000  -0.5000000  -0.5000000   0.0117188  crystal
kpoint   0.1250000   0.2500000   0.3750000   0.0937500  crystal
kpoint   0.1250000   0.2500000  -0.5000000   0.1875000  crystal
kpoint   0.1250000   0.2500000  -0.3750000   0.0937500  crystal
kpoint   0.1250000   0.3750000  -0.5000000   0.0937500  crystal
kpoint   0.1250000   0.3750000  -0.3750000   0.1875000  crystal
kpoint   0.1250000   0.3750000  -0.2500000   0.0937500  crystal
kpoint   0.1250000  -0.5000000  -0.3750000   0.0468750  crystal
kpoint   0.2500000  -0.5000000  -0.2500000   0.0234375  crystal

set threshold_scf 1.E-10 50

randomize_wf
load -states wfsave/si54p.wf
run 0 300 2
save -states wfsave/si54p.wfb
run 0 300 2
save -states wfsave/si54p.wfc
run 0 300 2
save -states wfsave/si54p.wfd

quit
