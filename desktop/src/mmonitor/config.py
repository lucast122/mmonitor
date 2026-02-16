from os.path import dirname, abspath, join

# STATIC PATHS  – resources now live inside the package
_MMONITOR_ROOT = dirname(abspath(__file__))
_RESOURCES = join(_MMONITOR_ROOT, 'resources')
images_path = join(_RESOURCES, 'images')
r_path = join(_RESOURCES, 'r')
