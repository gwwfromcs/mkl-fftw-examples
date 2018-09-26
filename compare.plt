
#set logscale y

plot '30_3_2_1e-5.dat' u ($1):(abs($6)) w lp lc rgb "#5dc19b",\
     '40_3_2_1e-5.dat' u ($1):(abs($6)) w lp lc rgb "#57c8cc",\
     '50_3_2_1e-5.dat' u ($1):(abs($6)) w lp lc rgb "#57accc",\
     '60_3_2_1e-5.dat' u ($1):(abs($6)) w lp lc rgb "#5789cc",\
     '70_3_2_1e-5.dat' u ($1):(abs($6)) w lp lc rgb "#576acc",\
     '80_3_2_1e-5.dat' u ($1):(abs($6)) w lp lc rgb "#7d57cc",\
    '120_3_2_1e-5.dat' u ($1):(abs($6)) w lp lc rgb "#b857cc",\
    '160_3_2_1e-5.dat' u ($1):(abs($6)) w lp lc rgb "#cc577d",\
    '200_3_2_1e-5.dat' u ($1):(abs($6)) w lp lc rgb "#cc7f57",\
    '300_3_2_1e-5.dat' u ($1):(abs($6)) w lp lc rgb "#95cc57",\
    '500_3_2_1e-5.dat' u ($1):(abs($6)) w lp lc rgb "#c2cc57"
