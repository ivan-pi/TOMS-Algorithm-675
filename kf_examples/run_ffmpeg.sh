# https://stackoverflow.com/questions/24961127/how-to-create-a-video-from-images-with-ffmpeg
ffmpeg -framerate 30 -i img%05d.png -c:v libx264 -pix_fmt yuv420p drying2.mp4
