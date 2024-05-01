import os


_file_dir = os.path.dirname(__file__)
_parent_dir = os.path.dirname(_file_dir)

nuke.pluginAddPath(_file_dir)
nuke.pluginAddPath(os.path.join(_parent_dir, "gizmos"))
