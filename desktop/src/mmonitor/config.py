from os.path import dirname, abspath, join

# STATIC PATHS
_MMONITOR_ROOT = dirname(abspath(__file__))
_RESOURCES = join(dirname(_MMONITOR_ROOT), 'resources')
images_path = join(_RESOURCES, 'images')
r_path = join(_RESOURCES, 'r')
