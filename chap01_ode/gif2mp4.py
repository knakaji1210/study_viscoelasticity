import sys
import moviepy.editor as mp

# 使い方
# python3 gif2mp4.py ./gif/Brownian_Motion_100steps.gif

args = sys.argv

if __name__ == '__main__':
    input_file_name = args[1]
    output_file_name = input_file_name[0:-3]+"mp4"
    clip = mp.VideoFileClip(input_file_name)
    clip.write_videofile(output_file_name)

