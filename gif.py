import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import os
import pandas as pd
from PIL import Image
import re



def create_gif(image_dir, output_path, image_suffix, image_step=1):
    images = []
    t = 0
    while True:
        image_file = os.path.join(image_dir, f"t{t}_{image_suffix}.png")
        if not os.path.exists(image_file):
            break

        print(f"Adding image {image_file}")
        img = Image.open(image_file)
        images.append(img)

        t += image_step
        if t > 3000:
            break

    if images:
        # Save the first image as the GIF background
        images[0].save(output_path, save_all=True, append_images=images[1:], loop=0, duration=image_step*10)
        print(f"GIF created successfully at {output_path}")
    else:
        print("No images found to create the GIF.")


local = False
repo = '20230705_lk001_Linux/CONFIG_vasculature_irrad_single_example.py_162736'

iter = [3]
image_step = 4


csv_file = ''
#all repositories in repo:
#find csv file in repo
for filename in os.listdir(repo):
    # Check if the file is a csv file
    if filename.endswith('.csv'):
        csv_file = filename

param_space = pd.read_csv(f'{repo}/{csv_file}', sep=' ', header=0)
number_of_iterations = len(param_space['Iteration'])

paths = [f'{repo}/iter{i}/DataOutput/' for i in range(0, number_of_iterations)]
#remove paths 4
print(paths)


if local: paths = ['DataOutput/']


for i in iter:
    print('iter', i, 'of', iter, 'is being processed')
    image_directory = repo + '/iter' + str(i) + '/Plots/CurrentPlotting'

    image_sufix1 = 'AllPlots'
    output_path1 = repo + '/' + image_sufix1 + str(i) + '.gif'
    image_sufix2 = 'Vasculature'
    output_path2 = repo + '/' + image_sufix2 + str(i) + '.gif'
    if local:
        image_directory = 'Plots/CurrentPlotting'
        output_path1 = image_sufix1 + '.gif'
        output_path2 = image_sufix2 + '.gif'

    create_gif(image_directory, output_path1, image_sufix1, image_step)
    create_gif(image_directory, output_path2, image_sufix2, image_step)


