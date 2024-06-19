# Blink Fog 3D

This gizmo allows you to create noise in 3D space using just a camera as input.

![fog_ex](https://github.com/obulka/blink_fog_3d/assets/21975584/aeeb19b2-c347-47d2-9d7d-919786fe870c)

blink_fog_3d is a tool for rendering 4d noise through a render camera with lots of additional user options to change the look and feel of the noise. The advantage of this blinkscript versus traditional 3d noise setups is the noise is rendered entirely in blinkscript and doesn't require the user to place a box or cards of noise in the scene. The node has a built in depth ramp so you can choose what depth from the camera the noise is visible and choose how the noise falls off as it gets closer to the camera (useful for making fast moving wispies). The 4d aspect of the noise allows you to animate the seed of the noise without having a visible direction to the noise (xyz) which can be used to make the noise look like smoke or steam. A depth input can also be used to apply holdouts to the noise for use as a procedural fog generator. Note that there is an 'invert depth' checkbox to use 1/z. The samples per ray knob is the general quality slider, the more samples you provide, the less noise you get but the slower the node becomes. The node is pretty fast if you have a fast gpu and a $GUI expression can be used on the samples per ray knob to make it more responsive locally.

## Usage

You can add the `src/python` directory to your `NUKE_PATH` to add the gizmo to your node menu or simply open the example file to find the group node.

You can also pass a depth AOV to the node to create holdouts.

### New in Version 1.4.0

- Added a `y ramp` knob to vary the noise density along the y-axis in worldspace

### New in Version 1.3.0

- Separate deep holdout input because nuke does not like mixing deep and 2d data

### New in Version 1.2.0

- Deep support!
  - Use the `holdout mode` `deep (proxy)` for fast deep that will work for most scenarios, but can have artifacts
  - Use the `holdout mode` `deep (full)` for slower but completely accurate deep
  - Use the `output deep` checkbox to output deep data without holdouts
 
- Overscan
  - Use the `overscan` knob to specify any overscan
