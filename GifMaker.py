# Makes a gif of desired images

import glob
import os
from PIL import Image
def make_gif(frame_folder):

    files = (glob.glob(f"{frame_folder}/*.jpg"))

    
    files_ascending = sorted(files, key=lambda t: os.stat(t).st_mtime)
    print(files_ascending)
    #files_descending = sorted(files, key=lambda t: -os.stat(t).st_mtime)
    #print(files_descending)
    frames = [Image.open(image) for image in files_ascending]
    
    frame_one = (frames)[0]
    frame_one.save("turtles.gif", format="GIF", append_images=frames,
               save_all=True, duration=1000, loop=0)
    
if __name__ == "__main__":
    make_gif("Images")