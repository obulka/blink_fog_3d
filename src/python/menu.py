import nuke

_toolbar = nuke.toolbar("Nodes")

_menu = _toolbar.addMenu("Blink3D")
_menu.addCommand("blink_fog_3d", "nuke.createNode('blink_fog_3d')")
