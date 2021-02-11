gata_counts = read.table('count.txt', header=T)
gata_counts$region = NULL

no_tgataa = gata_counts[ gata_counts[[1]] == 0,]
no_wgatar = no_tgataa[ no_tgataa[[2]] == 0,]

bp = barplot(table(gata_counts[[1]]), main="Breakdown of TGATAA-containing peaks")
text(bp, rep(250,12),labels=table(gata_counts[[1]]),srt=90)
bp = barplot(table(no_tgataa[[2]]), main="Breakdown of non-TGATAA, WGATAR-containing peaks")
text(bp, rep(100,10),labels=table(no_tgataa[[2]]),srt=90)
bp = barplot(table(no_wgatar[[3]]), main="Breakdown of non-WGATAR, GATA-containing peaks")
text(bp, rep(25,14),labels=table(no_wgatar[[3]]),srt=90)