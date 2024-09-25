# Blink Fog 3D

This gizmo allows you to create noise in 3D space using just a camera as input.

![fog_ex](https://github.com/obulka/blink_fog_3d/assets/21975584/aeeb19b2-c347-47d2-9d7d-919786fe870c)

blink_fog_3d is a tool for rendering 4d noise through a render camera with lots of additional user options to change the look and feel of the noise. The advantage of this blinkscript versus traditional 3d noise setups is the noise is rendered entirely in blinkscript and doesn't require the user to place a box or cards of noise in the scene. The node has a built in depth ramp so you can choose what depth from the camera the noise is visible and choose how the noise falls off as it gets closer to the camera (useful for making fast moving wispies). The 4d aspect of the noise allows you to animate the seed of the noise without having a visible direction to the noise (xyz) which can be used to make the noise look like smoke or steam. A depth input can also be used to apply holdouts to the noise for use as a procedural fog generator. Note that there is an 'invert depth' checkbox to use 1/z. The samples per ray knob is the general quality slider, the more samples you provide, the less noise you get but the slower the node becomes. The node is pretty fast if you have a fast gpu and a $GUI expression can be used on the samples per ray knob to make it more responsive locally.

## Usage

You can add the `src/python` directory to your `NUKE_PATH` to add the gizmo to your node menu or simply open the example file to find the group node.

You can also pass a depth AOV to the node to create holdouts.

### New in Version 2.3.0

- Removed the 'pixel offset' knob and inherit that functionality from the camera input
  - Added support for 'window translate' and 'window scale' knobs on the camera input

### New in Version 2.2.0

- Added a 'pixel offset' knob to effectively translate the image window before rendering, this can be animated to match camera shake.
- Fixed a bug where shapes would dim when depth ramp enabled because sampling started at beginning of depth range

### New in Version 2.1.0

- Added a 'max iterations' knob to further restrict the number of ray marching steps
- Fixed bug when objects were > 100000 units away from camera which, due to loss of precision, would cause the normal to be nan
- Ensured that 'hit tolerance' is positive

### New in Version 2.0.0

- Depth ramp is now optional
  - This allows you to use the shape options fully in world space without worrying about them leaving the depth range during a camera move
  - If shapes are used we now start sampling noise where the ray first intersects the shape and stop where the ray exits the shape
    - This enables much more optimal sampling
- The various ramp groups (other than depth) have been simplified into a single dropdown menu
- The planar ramp has been replace with 'in', 'out', and 'falloff' knobs
  - All of these values are expected to be a positive distance from the plane, where in is the distance below, and out is the distance above
  - The distance from the plane can now be specified on either side individually rather than being mirrored
- If the depth was in camera space and the holdout was far from the camera and moved from the center of frame to the edge, it would appear to move in front of the fog
  - Added a 'camera space' checkbox so that if the depth is in camera space the sampling can account for this
- Ensured that node and example work in Nuke 12 and 15

### New in Version 1.8.0

- Added exponential falloff to sphere, plane, and box ramps

### New in Version 1.7.0

- Added a `box ramp` to contain the noise within a box at a position and with a rotation

### New in Version 1.6.0

- Added a `spherical ramp` allowing you to create a sphere of noise in the scene with a linear falloff
  - If this and the `planar ramp` are both enabled, the intersection of the two will be used

### New in Version 1.5.0

- Changed the `y ramp` to a `planar ramp`
  - allows you to specify a plane of noise with a falloff on either side
  - the plane is specified by a point and normal, allowing any arbitrary direction/positioning of the noise plane

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
