import os
import shutil
from os.path import join, dirname, realpath
from subprocess import call

ROOT = dirname(realpath(__file__))
# if getattr(sys, 'frozen', False):
#     # If it's bundled, adjust the root path to the directory where the executable resides
#     ROOT = sys._MEIPASS
# else:
#     ROOT = os.path.dirname(os.path.realpath(__file__))

# Base path for images
IMAGES_PATH = os.path.join(ROOT, "src", "resources", "images")



def main():
    entry_point = join(ROOT, 'src', 'mmonitor', '__main__.py')
    placeholder = join(ROOT, 'src', 'resources', 'images', '.placeholder')
    images_dir_dest = join('resources', 'images')
    r_script = join(ROOT, 'src', 'resources', 'r', 'horizon.r')
    r_dir_dest = join('resources', 'r')

    # call(f"pyinstaller -D {entry_point} --name mmonitor --add-data {placeholder}{pathsep}{images_dir_dest} --add-data {r_script}{pathsep}{r_dir_dest}".split()) #use this for building new spec file
    # remove old dist and built for faster rebuilding
    try:
        shutil.rmtree("/Users/timolucas/PycharmProjects/MMonitor/desktop/build/mmonitor/")
        shutil.rmtree("/Users/timolucas/PycharmProjects/MMonitor/desktop/dist/mmonitor/")
    except:
        call(f"pyinstaller mmonitor.spec".split())  # use this for using edited spec file
    call(f"pyinstaller mmonitor.spec".split())  # use this for using edited spec file


if __name__ == '__main__':
    main()