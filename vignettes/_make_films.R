dirs <- list.dirs()
dir_bak <- getwd()

cmd_base <- "ffmpeg -f image2 -r 2 -i %d.jpg -vcodec libx264 -profile:v high444 -crf 0 "
for (dir in dirs){
	setwd(dir)
	print(dir)
	cmd <- paste0(cmd_base, dir, ".mp4")
	system(cmd)
	setwd(dir_bak)
}

