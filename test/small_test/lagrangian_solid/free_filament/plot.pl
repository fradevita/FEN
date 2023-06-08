reset
set xr [-0.4:0.4]
set yr [-1:0]
set size ratio 1
do for [ii = 0:32000:500] {
    filename = sprintf('data/sb_%07d', ii)
    plot 'data/sb_0000000' w l, filename w l t sprintf('t = %f', ii/100)
    pause 0.5
}
