import pathlib
cwd = pathlib.Path.cwd()
import TimeEmbedding as TE
from PIL import Image, ImageFont, ImageDraw
import os
import shutil



def main():
    '''
    '''

    k_p = 9
    k_m = 9
    max_rows = 1000
    gap_min, gap_max = 0, 20

    image_dir = cwd / "TE_Images"
    exon_dir = image_dir / "Exon"
    intron_dir = image_dir / "Intron"
    both_dir = image_dir / "Both"
    
    for gap in range(gap_min, gap_max):
        print(f"### Gap = {gap} ###")
        TE.time_embedding(k_p = k_p, k_m = k_m, max_rows = max_rows, gap = gap)

    generator(exon_dir, str(cwd / f"Exons_gap_{gap_min}_{gap_max}.gif"))
    generator(intron_dir, str(cwd / f"Introns_gap_{gap_min}_{gap_max}.gif"))
    generator(both_dir, str(cwd / f"Both_gap_{gap_min}_{gap_max}.gif"))


def generator(specific_folder: pathlib.Path, output_file):
    '''
    '''
    temp_folder: pathlib.Path = cwd / "TempPNGImages"
    if not temp_folder.is_dir():
        temp_folder.mkdir()

    png_files = list(specific_folder.rglob("*.png"))
    png_files.sort(key = lambda x: os.path.getmtime(x))

    frames = []
    for png in png_files:
        new_frame = Image.open(png)
        # print(png)
        frames.append(new_frame)


    frames[0].save(output_file, format="GIF", append_images=frames[1:], save_all = True, loop = 0, duration = 225)

    shutil.rmtree(str(temp_folder))



if __name__ in "__main__":
    main()