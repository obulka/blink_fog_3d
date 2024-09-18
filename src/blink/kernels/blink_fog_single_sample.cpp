// Copyright 2024 by Owen Bulka.
// All rights reserved.
// This file is released under the "MIT License Agreement".
// Please see the LICENSE.md file that should have been included as part
// of this package.


//
// Math operations
//


/**
 * The maximum component of a vector.
 *
 * @arg vector: The vector.
 *
 * @returns: The maximum component of the vector.
 */
inline float maxComponent(const float3 &vector)
{
    return max(vector.x, max(vector.y, vector.z));
}


/**
 * The negative part of the vector. Ie. any positive values will be 0,
 * and the negative values will be positive.
 *
 * @arg vector: The vector.
 *
 * @returns: The negative part of the vector.
 */
inline float negativePart(const float value)
{
    return -min(value, 0.0f);
}


/**
 * The positive part of the vector. Ie. any negative values will be 0.
 *
 * @arg vector: The vector.
 *
 * @returns: The positive part of the vector.
 */
inline float3 positivePart(const float3 &vector)
{
    return max(vector, float3(0));
}


/**
 * Compute the fractional portion of the value. Ex. 3.5 returns 0.5
 *
 * @arg value: The value to get the fractional portion of.
 *
 * @returns: The fractional portion of the value.
 */
inline float fract(const float value)
{
    return value - floor(value);
}


/**
 * Saturate a value ie. clamp between 0 and 1
 *
 * @arg value: The value to saturate
 *
 * @returns: The clamped value
 */
inline float saturate(float value)
{
    return clamp(value, 0.0f, 1.0f);
}


/**
 * Multiply a 4d vector by a 4x4 matrix.
 *
 * @arg m: The matrix that will transform the vector.
 * @arg v: The vector to transform.
 * @arg out: The location to store the resulting vector.
 */
inline void matmul(const float4x4 &m, const float4 &v, float4 &out)
{
    for (int i=0; i < 4; i++)
    {
        out[i] = 0;

        for (int j=0; j < 4; j++)
        {
            out[i] += m[i][j] * v[j];
        }
    }
}


/**
 * Multiply a 4d vector by a 4x4 matrix.
 *
 * @arg m: The matrix that will transform the vector.
 * @arg v: The vector to transform.
 * @arg out: The location to store the resulting vector.
 */
inline float4 matmul(const float4x4 &m, const float4 &v)
{
    float4 out;
    for (int i=0; i < 4; i++)
    {
        out[i] = 0;

        for (int j=0; j < 4; j++)
        {
            out[i] += m[i][j] * v[j];
        }
    }

    return out;
}


/**
 * Multiply a 3d vector by a 3x3 matrix.
 *
 * @arg m: The matrix that will transform the vector.
 * @arg v: The vector to transform.
 * @arg out: The location to store the resulting vector.
 */
inline float3 matmul(const float3x3 &m, const float3 &v)
{
    return float3(
        m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
        m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
        m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z
    );
}


/**
 * Convert the uv position in a latlong image to angles.
 *
 * @arg uvPosition: The UV position.
 *
 * @returns: The equivalent angles in radians.
 */
inline float2 uvPositionToAngles(const float2 &uvPosition)
{
    return float2(
        (uvPosition.x + 1.0f) * PI,
        (1.0f - uvPosition.y) * PI / 2.0f
    );
}


/**
 * Convert location of a pixel in an image into UV.
 *
 * @arg pixelLocation: The x, and y positions of the pixel.
 * @arg format: The image width, and height.
 *
 * @returns: The UV position.
 */
inline float2 pixelsToUV(const float2 &pixelLocation, const float2 &format)
{
    return float2(
        2.0f * pixelLocation.x / format.x - 1.0f,
        2.0f * pixelLocation.y / format.y - 1.0f
    );
}


/**
 * Convert a spherical unit vector (unit radius) to cartesion.
 *
 * @arg angles: The spherical angles in radians.
 *
 * @returns: The equivalent cartesion vector.
 */
inline float3 sphericalUnitVectorToCartesion(const float2 &angles)
{
    const float sinPhi = sin(angles.y);
    return float3(
        cos(angles.x) * sinPhi,
        cos(angles.y),
        sin(angles.x) * sinPhi
    );
}


/**
 * Get the position component of a world matrix.
 *
 * @arg worldMatrix: The world matrix.
 * @arg position: The location to store the position.
 */
inline void positionFromWorldMatrix(const float4x4 &worldMatrix, float3 &position)
{
    position = float3(
        worldMatrix[0][3],
        worldMatrix[1][3],
        worldMatrix[2][3]
    );
}


/**
 * Get the rotation component of a world matrix.
 *
 * @arg worldMatrix: The world matrix.
 * @arg rotation: The location to store the rotation.
 */
inline void rotationFromWorldMatrix(const float4x4 &worldMatrix, float3x3 &rotationMatrix)
{
    rotationMatrix[0][0] = worldMatrix[0][0];
    rotationMatrix[0][1] = worldMatrix[0][1];
    rotationMatrix[0][2] = worldMatrix[0][2];
    rotationMatrix[1][0] = worldMatrix[1][0];
    rotationMatrix[1][1] = worldMatrix[1][1];
    rotationMatrix[1][2] = worldMatrix[1][2];
    rotationMatrix[2][0] = worldMatrix[2][0];
    rotationMatrix[2][1] = worldMatrix[2][1];
    rotationMatrix[2][2] = worldMatrix[2][2];
}


/**
 * Get a rotation matrix from radian angle values in ZYX order.
 *
 * @arg angles: The rotation angles in radians.
 * @arg out: The location to store the rotation matrix.
 */
inline void rotationMatrix(const float3 &angles, float3x3 &out)
{
    const float3 cosAngles = cos(angles);
    const float3 sinAngles = sin(angles);

    // Why tf can I not init a float3x3 normally??
    out[0][0] = cosAngles.y * cosAngles.z;
    out[0][1] = sinAngles.x * sinAngles.y * cosAngles.z - cosAngles.x * sinAngles.z;
    out[0][2] = cosAngles.x * sinAngles.y * cosAngles.z + sinAngles.x * sinAngles.z;
    out[1][0] = cosAngles.y * sinAngles.z;
    out[1][1] = sinAngles.x * sinAngles.y * sinAngles.z + cosAngles.x * cosAngles.z;
    out[1][2] = cosAngles.x * sinAngles.y * sinAngles.z - sinAngles.x * cosAngles.z;
    out[2][0] = -sinAngles.y;
    out[2][1] = sinAngles.x * cosAngles.y;
    out[2][2] = cosAngles.x * cosAngles.y;
}


//
// Randomization functions
//


// Some random constants on the interval [1, 2]
#define RAND_CONST_0 1.571411510193971f


/**
 * Get a random value on the interval [0, 1].
 *
 * @arg seed: The random seed.
 *
 * @returns: A random value on the interval [0, 1].
 */
inline float random(const float seed)
{
    return fract(sin(73.1f * seed + 91.3458f) * 47453.5453f);
}


/**
 * Get a random value on the interval [0, 1].
 *
 * @arg seed: The random seed.
 *
 * @returns: A random value on the interval [0, 1].
 */
inline float2 random(const float2 &seed)
{
    return float2(
        fract(sin(13.157f * seed.x + 71.743f) * 7513.471f),
        fract(sin(97.519f * seed.y + 113.591f) * 47453.5453f)
    );
}


/**
 * Create a random point that lies within the unit circle.
 *
 * @arg seed: The random seed.
 *
 * @returns: A random point, (radius, angle) in the unit circle.
 */
inline float2 uniformPointInUnitCircle(const float2 &seed)
{
    return float2(sqrt(random(seed.x)), 2.0f * PI * random(seed.y));
}


/**
 * Sum the elements of a vector.
 *
 * @arg vector: The vector to sum the elements of.
 *
 * @returns: The sum of the vector elements.
 */
inline float elementSum(const float4 &vector) {
    return vector.x + vector.y + vector.z + vector.w;
}


#define FBM_NOISE 1


//
// Camera Utilities
//


/**
 * Compute the field of view from focal length.
 *
 * @arg focalLength: The focal length.
 *
 * @returns: The equivalent field of view.
 */
inline float fieldOfView(const float focalLength)
{
    return 2 * atan(1 / focalLength);
}


/**
 * Compute the aspect ratio from image format.
 *
 * @arg height_: The height of the image.
 * @arg width_: The width of the image.
 *
 * @returns: The aspect ratio.
 */
inline float aspectRatio(const float height_, const float width_)
{
    return height_ / width_;
}


/**
 *
 */
inline float fStopToAperture(const float fStop, const float focalLength)
{
    return focalLength / fStop / 1000.0f;
}


/**
 * Create a projection matrix for a camera.
 *
 * @arg focalLength: The focal length of the camera.
 * @arg horizontalAperture: The horizontal aperture of the camera.
 * @arg aspect: The aspect ratio of the camera.
 * @arg nearPlane: The distance to the near plane of the camera.
 * @arg farPlane: The distance to the far plane of the camera.
 *
 * @returns: The camera's projection matrix.
 */
float4x4 projectionMatrix(
        const float focalLength,
        const float horizontalAperture,
        const float aspect,
        const float nearPlane,
        const float farPlane)
{
    float farMinusNear = farPlane - nearPlane;
    return float4x4(
        2 * focalLength / horizontalAperture, 0, 0, 0,
        0, 2 * focalLength / horizontalAperture / aspect, 0, 0,
        0, 0, -(farPlane + nearPlane) / farMinusNear, -2 * (farPlane * nearPlane) / farMinusNear,
        0, 0, -1, 0
    );
}


/**
 * Generate a ray out of a camera.
 *
 * @arg cameraWorldMatrix: The camera matrix.
 * @arg inverseProjectionMatrix: The inverse of the projection matrix.
 * @arg uvPosition: The UV position in the resulting image.
 * @arg rayOrigin: Will store the origin of the ray.
 * @arg rayDirection: Will store the direction of the ray.
 */
void createCameraRay(
        const float4x4 &cameraWorldMatrix,
        const float4x4 &inverseProjectionMatrix,
        const float2 &uvPosition,
        float3 &rayOrigin,
        float3 &rayDirection)
{
    positionFromWorldMatrix(cameraWorldMatrix, rayOrigin);
    float4 direction = matmul(
        inverseProjectionMatrix,
        float4(uvPosition.x, uvPosition.y, 0, 1)
    );
    matmul(
        cameraWorldMatrix,
        float4(direction.x, direction.y, direction.z, 0),
        direction
    );
    rayDirection = normalize(float3(direction.x, direction.y, direction.z));
}


/**
 * Generate a ray out of a camera.
 *
 * @arg cameraWorldMatrix: The camera matrix.
 * @arg inverseProjectionMatrix: The inverse of the projection matrix.
 * @arg uvPosition: The UV position in the resulting image.
 * @arg rayOrigin: Will store the origin of the ray.
 * @arg rayDirection: Will store the direction of the ray.
 */
void createCameraRay(
        const float4x4 &cameraWorldMatrix,
        const float4x4 &inverseProjectionMatrix,
        const float2 &uvPosition,
        const float aperture,
        const float focalDistance,
        const float2 &seed,
        float3 &rayOrigin,
        float3 &rayDirection)
{
    createCameraRay(
        cameraWorldMatrix,
        inverseProjectionMatrix,
        uvPosition,
        rayOrigin,
        rayDirection
    );

    const float4 cameraForward4 = matmul(
        cameraWorldMatrix,
        float4(0, 0, -1, 0)
    );
    const float3 cameraForward = float3(
        cameraForward4.x,
        cameraForward4.y,
        cameraForward4.z
    );
    const float4 cameraRight4 = matmul(
        cameraWorldMatrix,
        float4(1, 0, 0, 0)
    );
    const float3 cameraRight = float3(
        cameraRight4.x,
        cameraRight4.y,
        cameraRight4.z
    );
    const float4 cameraUp4 = matmul(
        cameraWorldMatrix,
        float4(0, 1, 0, 0)
    );
    const float3 cameraUp = float3(
        cameraUp4.x,
        cameraUp4.y,
        cameraUp4.z
    );

    const float3 focalPlanePoint = rayOrigin + cameraForward * focalDistance;
    const float3 focalPlaneNormal = -cameraForward;

    const float focalPointDistance = (
        (dot(focalPlaneNormal, focalPlanePoint) - dot(rayOrigin, focalPlaneNormal))
        / dot(rayDirection, focalPlaneNormal)
    );
    const float3 focalPoint = rayOrigin + focalPointDistance * rayDirection;

    const float2 pointInUnitCircle = uniformPointInUnitCircle(seed);
    const float2 offset = pointInUnitCircle.x * aperture * float2(
        cos(pointInUnitCircle.y),
        sin(pointInUnitCircle.y)
    );

    rayOrigin += cameraRight * offset.x + cameraUp * offset.y;
    rayDirection = normalize(focalPoint - rayOrigin);
}


/**
 * Generate a LatLong ray out of a camera.
 *
 * @arg cameraWorldMatrix: The camera matrix.
 * @arg uvPosition: The UV position in the resulting image.
 * @arg rayOrigin: Will store the origin of the ray.
 * @arg rayDirection: Will store the direction of the ray.
 */
void createLatLongCameraRay(
        const float4x4 &cameraWorldMatrix,
        const float2 &uvPosition,
        float3 &rayOrigin,
        float3 &rayDirection)
{
    positionFromWorldMatrix(cameraWorldMatrix, rayOrigin);
    rayDirection = sphericalUnitVectorToCartesion(uvPositionToAngles(uvPosition));

    float3x3 cameraRotation;
    rotationFromWorldMatrix(cameraWorldMatrix, cameraRotation);
    rayDirection = matmul(cameraRotation, rayDirection);
}


// SDFs


/**
 * Compute the signed distance along a vector
 *
 * @arg vector: A vector from a point to the nearest surface of an
 *     object.
 *
 * @returns: The signed length of the vector.
 */
inline float sdfLength(const float3 &vector)
{
    return (
        length(positivePart(vector))
        - negativePart(maxComponent(vector))
    );
}


/**
 * Compute the min distance from a point to a sphere.
 *
 * @arg position: The point to get the distance to, from the object.
 * @arg radius: The radius of the sphere.
 *
 * @returns: The minimum distance from the point to the shape.
 */
inline float distanceToSphere(const float3 &position, const float radius)
{
    return length(position) - radius;
}


/**
 * Compute the min distance from a point to a plane.
 * Anything underneath the plane, as defined by the normal direction
 * pointing above, will be considered inside.
 *
 * @arg position: The point to get the distance to, from the object.
 * @arg normal: The normal direction of the plane.
 *
 * @returns: The minimum distance from the point to the shape.
 */
inline float distanceToPlane(const float3 &position, const float3 &normal)
{
    return dot(position, normal);
}


/**
 * Compute the min distance from a point to a rectangular prism.
 * Centered at the origin.
 *
 * @arg position: The point to get the distance to, from the object.
 * @arg width: The width (x) of the prism.
 * @arg height: The height (y) of the prism.
 * @arg depth: The depth (z) of the prism.
 *
 * @returns: The minimum distance from the point to the shape.
 */
inline float distanceToRectangularPrism(
        const float3 &position,
        const float width,
        const float height,
        const float depth)
{
    // Only look at positive quadrant, using symmetry
    const float3 prismToPosition = fabs(position) - float3(width, height, depth) / 2.0f;
    // Clamp the components that are inside the prism to the surface
    // before getting the distance
    return sdfLength(prismToPosition);
}


kernel SingleSampleFogKernel : ImageComputationKernel<ePixelWise>
{
    // the input which specifies the format, process is called once per pixel
    // in this image, which also provides random seeds
    Image<eRead, eAccessPoint, eEdgeNone> seeds;

    // the output image
    Image<eWrite> dst;


    param:
        // These parameters are made available to the user.

        // Camera params
        float _focalLength;
        float _horizontalAperture;
        float _nearPlane;
        float _farPlane;
        float4x4 _cameraWorldMatrix;
        float _focalDistance;
        float _fStop;
        bool _depthOfFieldEnabled;
        bool _latLong;
        bool _useCameraDepth;

        // Image params
        float _formatWidth;
        float _formatHeight;
        float _individualSample;
        float _density;
        int _samplesPerRay;

        bool _enableDepthRamp;
        float4 _depthRamp;

        float4 _planarRamp;
        bool _enablePlanarRamp;
        float3 _planePosition;
        float3 _planeNormal;
        float _planarFalloffPower;
        float _planarFalloffOffset;

        bool _enableSphericalRamp;
        float3 _sphericalRampPosition;
        float2 _sphericalRampRadii;
        float _sphericalFalloffPower;
        float _sphericalFalloffOffset;

        float3 _boxRampPosition;
        float3 _boxRampRotation;
        float3 _boxRampDimensions;
        float _boxRampFalloff;
        bool _enableBoxRamp;
        float _boxFalloffPower;
        float _boxFalloffOffset;

        bool _secondSample;

        // Noise Parameters
        float _size;
        int _noiseType;
        float3 _translation;
        float _octaves;
        float _lacunarity;
        float _gain;
        float _gamma;
        float4 _lowFrequencyScale;
        float4 _highFrequencyScale;
        float4 _lowFrequencyTranslation;
        float4 _highFrequencyTranslation;

        // Raymarch parameters
        float _hitTolerance;
        float _maxDistance;

    local:
        // These local variables are not exposed to the user.
        float4x4 __inverseCameraProjectionMatrix;
        float4x4 __inverseCameraWorldMatrix;
        float __aperture;

        float4 __translation;

        int __simplex[64][4];
        int __perm[512];
        int __grad4[32][4];

        float __depthSampleStep;

        float __sphericalRampFalloffDistance;

        float3 __planeNormal;

        float3x3 __inverseBoxRotMatrix;

        bool __enableRaymarching;

        float3 __offset0;
        float3 __offset1;
        float3 __offset2;
        float3 __offset3;


    /**
     * Give the parameters labels and default values.
     */
    void define()
    {
        // Camera params
        defineParam(_focalLength, "FocalLength", 50.0f);
        defineParam(_horizontalAperture, "HorizontalAperture", 24.576f);
        defineParam(_nearPlane, "NearPlane", 0.1f);
        defineParam(_farPlane, "FarPlane", 10000.0f);
        defineParam(
            _cameraWorldMatrix,
            "CameraWorldMatrix",
            float4x4(
                1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1
            )
        );
        defineParam(_focalDistance, "FocalDistance", 4.0f);
        defineParam(_fStop, "fstop", 16.0f);
        defineParam(_depthOfFieldEnabled, "EnableDepthOfField", true);
        defineParam(_latLong, "OutputLatLong", false);
        defineParam(_useCameraDepth, "UseCameraDepth", false);

        // Image params
        defineParam(_formatHeight, "ScreenHeight", 2160.0f);
        defineParam(_formatWidth, "ScreenWidth", 3840.0f);
        defineParam(_individualSample, "IndividualSample", 0.0f);
        defineParam(_density, "Density", 1.0f);
        defineParam(_samplesPerRay, "SamplesPerRay", 5);

        defineParam(_enableDepthRamp, "EnableDepthRamp", true);
        defineParam(_depthRamp, "SampleDepthRamp", float4(5.0f, 10.0f, 15.0f, 20.0f));

        defineParam(_planarRamp, "SamplePlanarRamp", float4(5.0f, 10.0f, 15.0f, 20.0f));
        defineParam(_enablePlanarRamp, "EnablePlanarRamp", false);
        defineParam(_planePosition, "PlanePosition", float3(0.0f, 0.0f, 0.0f));
        defineParam(_planeNormal, "PlaneNormal", float3(0.0f, 1.0f, 0.0f));
        defineParam(_planarFalloffPower, "PlanarFalloffPower", 0.0f);
        defineParam(_planarFalloffOffset, "PlanarFalloffOffset", 0.0f);

        defineParam(_enableSphericalRamp, "EnableSphericalRamp", false);
        defineParam(_sphericalRampPosition, "SphericalRampPosition", float3(0.0f, 0.0f, 0.0f));
        defineParam(_sphericalRampRadii, "SphericalRampRadii", float2(1.0f, 2.0f));
        defineParam(_sphericalFalloffPower, "SphericalFalloffPower", 0.0f);
        defineParam(_sphericalFalloffOffset, "SphericalFalloffOffset", 0.0f);

        defineParam(_boxRampPosition, "BoxRampPosition", float3(0.0f, 0.0f, 0.0f));
        defineParam(_boxRampRotation, "BoxRampRotation", float3(0.0f, 0.0f, 0.0f));
        defineParam(_boxRampDimensions, "BoxRampDimensions", float3(1.0f, 1.0f, 1.0f));
        defineParam(_boxRampFalloff, "BoxRampFalloff", 1.0f);
        defineParam(_enableBoxRamp, "EnableBoxRamp", false);
        defineParam(_boxFalloffPower, "BoxFalloffPower", 0.0f);
        defineParam(_boxFalloffOffset, "BoxFalloffOffset", 0.0f);

        defineParam(_secondSample, "DoSecondSample", false);

        defineParam(_hitTolerance, "HitTolerance", 0.001f);
        defineParam(_maxDistance, "MaxDistance", 1000000.0f);

        // Noise Parameters
        defineParam(_size, "Size", 20.0f);
        defineParam(_noiseType, "NoiseType", 0);
        defineParam(_translation, "Translation", float3(0.0f, 0.0f, 0.0f));
        defineParam(_octaves, "Octaves", 8.0f);
        defineParam(_lacunarity, "Lacunarity", 3.0f);
        defineParam(_gain, "Gain", 0.5f);
        defineParam(_gamma, "Gamma", 0.5f);
        defineParam(_lowFrequencyScale, "LowFrequencyScale", float4(0.0f, 0.0f, 0.0f, 0.0f));
        defineParam(_highFrequencyScale, "HighFrequencyScale", float4(0.0f, 0.0f, 0.0f, 0.0f));
        defineParam(_lowFrequencyTranslation, "LowFrequencyTranslation", float4(0.0f, 0.0f, 0.0f, 0.0f));
        defineParam(_highFrequencyTranslation, "HighFrequencyTranslation", float4(0.0f, 0.0f, 0.0f, 0.0f));
    }


    /**
     * Initialize the local variables.
     */
    void init()
    {
        const float aspect = aspectRatio(_formatHeight, _formatWidth);
        __inverseCameraProjectionMatrix = projectionMatrix(
            _focalLength,
            _horizontalAperture,
            aspect,
            _nearPlane,
            _farPlane
        ).invert();
        // Dumb af but trust that changing the following will either break nuke 12 or 15
        __inverseCameraWorldMatrix = _cameraWorldMatrix;
        __inverseCameraWorldMatrix = __inverseCameraWorldMatrix.invert();

        __aperture = fStopToAperture(_fStop, _focalLength);

        const int simplexInit[64][4] = {
            {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 0, 0, 0}, {0, 2, 3, 1},
            {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 2, 3, 0},
            {0, 2, 1, 3}, {0, 0, 0, 0}, {0, 3, 1, 2}, {0, 3, 2, 1},
            {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 3, 2, 0},
            {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
            {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
            {1, 2, 0, 3}, {0, 0, 0, 0}, {1, 3, 0, 2}, {0, 0, 0, 0},
            {0, 0, 0, 0}, {0, 0, 0, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
            {1, 0, 2, 3}, {1, 0, 3, 2}, {0, 0, 0, 0}, {0, 0, 0, 0},
            {0, 0, 0, 0}, {2, 0, 3, 1}, {0, 0, 0, 0}, {2, 1, 3, 0},
            {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
            {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
            {2, 0, 1, 3}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
            {3, 0, 1, 2}, {3, 0, 2, 1}, {0, 0, 0, 0}, {3, 1, 2, 0},
            {2, 1, 0, 3}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
            {3, 1, 0, 2}, {0, 0, 0, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}
        };

        for (int i = 0; i < 64; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                __simplex[i][j] = simplexInit[i][j];
            }
        }

        const int permInit[256] = {
            151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233,
            7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23,
            190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219,
            203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174,
            20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48,
            27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133,
            230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65,
            25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18,
            169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198,
            173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147,
            118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17,
            182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44,
            154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22,
            39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112,
            104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12,
            191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49,
            192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121,
            50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29,
            24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
        };

        for (int i = 0; i < 512; i++)
        {
            __perm[i] = permInit[i % 256];
        }

        const int grad4Init[32][4]= {
            {0, 1, 1, 1},  {0, 1, 1, -1},  {0, 1, -1, 1},  {0, 1, -1, -1},
            {0, -1, 1, 1}, {0, -1, 1, -1}, {0, -1, -1, 1}, {0, -1, -1, -1},
            {1, 0, 1, 1},  {1, 0, 1, -1},  {1, 0, -1, 1},  {1, 0, -1, -1},
            {-1, 0, 1, 1}, {-1, 0, 1, -1}, {-1, 0, -1, 1}, {-1, 0, -1, -1},
            {1, 1, 0, 1},  {1, 1, 0, -1},  {1, -1, 0, 1},  {1, -1, 0, -1},
            {-1, 1, 0, 1}, {-1, 1, 0, -1}, {-1, -1, 0, 1}, {-1, -1, 0, -1},
            {1, 1, 1, 0},  {1, 1, -1, 0},  {1, -1, 1, 0},  {1, -1, -1, 0},
            {-1, 1, 1, 0}, {-1, 1, -1, 0}, {-1, -1, 1, 0}, {-1, -1, -1, 0}
        };

        for (int i = 0; i < 32; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                __grad4[i][j] = grad4Init[i][j];
            }
        }

        __offset0 = 0.5773f * float3(1, -1, -1);
        __offset1 = 0.5773f * float3(-1, -1, 1);
        __offset2 = 0.5773f * float3(-1, 1, -1);
        __offset3 = 0.5773f * float3(1, 1, 1);

        __translation = float4(_translation.x, _translation.y, _translation.z, 0.0f);

        __depthSampleStep = (_depthRamp.w - _depthRamp.x) / (float) _samplesPerRay;

        __planeNormal = normalize(_planeNormal);

        __sphericalRampFalloffDistance = _sphericalRampRadii.y - _sphericalRampRadii.x;

        __enableRaymarching = _enableSphericalRamp || _enableBoxRamp || _enablePlanarRamp;

        rotationMatrix(_boxRampRotation, __inverseBoxRotMatrix);
        __inverseBoxRotMatrix.invert();
    }

    /**
     * 4D Perlin simplex noise
     *
     * Copyright (c) 2007-2012 Eliot Eshelman
     *
     * This program is free software: you can redistribute it and/or modify
     * it under the terms of the GNU General Public License as published by
     * the Free Software Foundation, either version 3 of the License, or
     * (at your option) any later version.
     *
     * This program is distributed in the hope that it will be useful,
     * but WITHOUT ANY WARRANTY; without even the implied warranty of
     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
     * GNU General Public License for more details.
     *
     * You should have received a copy of the GNU General Public License
     * along with this program. If not, see <http://www.gnu.org/licenses/>.
     *
     * @arg seed: The seed for the noise.
     *
     * @returns: Noise value in the range [-1, 1], value of 0 on all integer
     *     coordinates.
     */
    inline float perlinSimplexNoise(const float4 &seed) {
        const float G4 = 0.138196601f;
        const float4 i = floor(seed + elementSum(seed) * 0.309016994f);
        const float4 x0 = seed - i + elementSum(i) * G4;

        const int c = (
            ((x0.x > x0.y) << 5)
            | ((x0.x > x0.z) << 4)
            | ((x0.y > x0.z) << 3)
            | ((x0.x > x0.w) << 2)
            | ((x0.y > x0.w) << 1)
            | (x0.z > x0.w)
        );
        const int4 i1 = int4(
            __simplex[c][0] >= 3,
            __simplex[c][1] >= 3,
            __simplex[c][2] >= 3,
            __simplex[c][3] >= 3
        );
        const int4 i2 = int4(
            __simplex[c][0] >= 2,
            __simplex[c][1] >= 2,
            __simplex[c][2] >= 2,
            __simplex[c][3] >= 2
        );
        const int4 i3 = int4(
            __simplex[c][0] >= 1,
            __simplex[c][1] >= 1,
            __simplex[c][2] >= 1,
            __simplex[c][3] >= 1
        );

        const float4 x1 = x0 - float4(i1.x, i1.y, i1.z, i1.w) + G4;
        const float4 x2 = x0 - float4(i2.x, i2.y, i2.z, i2.w) + 2.0f * G4;
        const float4 x3 = x0 - float4(i3.x, i3.y, i3.z, i3.w) + 3.0f * G4;
        const float4 x4 = x0 - 1.0f + 4.0f * G4;

        const int ii = (int) i.x & 255;
        const int jj = (int) i.y & 255;
        const int kk = (int) i.z & 255;
        const int ll = (int) i.w & 255;

        const int gi0 = __perm[
            ii + __perm[jj + __perm[kk + __perm[ll]]]
        ] % 32;
        const int gi1 = __perm[
            ii + i1.x + __perm[jj + i1.y + __perm[kk + i1.z + __perm[ll + i1.w]]]
        ] % 32;
        const int gi2 = __perm[
            ii + i2.x + __perm[jj + i2.y + __perm[kk + i2.z + __perm[ll + i2.w]]]
        ] % 32;
        const int gi3 = __perm[
            ii + i3.x + __perm[jj + i3.y + __perm[kk + i3.z + __perm[ll + i3.w]]]
        ] % 32;
        const int gi4 = __perm[
            ii + 1 + __perm[jj + 1 + __perm[kk + 1 + __perm[ll + 1]]]
        ] % 32;

        float n0, n1, n2, n3, n4;
        float t0 = 0.6f - dot(x0, x0);
        if (t0 < 0)
        {
            n0 = 0.0f;
        }
        else
        {
            t0 *= t0;
            n0 = t0 * t0 * dot(
                float4(__grad4[gi0][0], __grad4[gi0][1], __grad4[gi0][2], __grad4[gi0][3]),
                x0
            );
        }

        float t1 = 0.6f - dot(x1, x1);
        if (t1 < 0)
        {
            n1 = 0.0f;
        }
        else
        {
            t1 *= t1;
            n1 = t1 * t1 * dot(
                float4(__grad4[gi1][0], __grad4[gi1][1], __grad4[gi1][2], __grad4[gi1][3]),
                x1
            );
        }

        float t2 = 0.6f - dot(x2, x2);
        if (t2 < 0)
        {
            n2 = 0.0f;
        }
        else
        {
            t2 *= t2;
            n2 = t2 * t2 * dot(
                float4(__grad4[gi2][0], __grad4[gi2][1], __grad4[gi2][2], __grad4[gi2][3]),
                x2
            );
        }

        float t3 = 0.6f - dot(x3, x3);
        if (t3 < 0)
        {
            n3 = 0.0f;
        }
        else
        {
            t3 *= t3;
            n3 = t3 * t3 * dot(
                float4(__grad4[gi3][0], __grad4[gi3][1], __grad4[gi3][2], __grad4[gi3][3]),
                x3
            );
        }

        float t4 = 0.6f - dot(x4, x4);
        if (t4 < 0)
        {
            n4 = 0.0f;
        }
        else {
            t4 *= t4;
            n4 = t4 * t4 * dot(
                float4(__grad4[gi4][0], __grad4[gi4][1], __grad4[gi4][2], __grad4[gi4][3]),
                x4
            );
        }

        return 27.0f * (n0 + n1 + n2 + n3 + n4);
    }


    //
    // Use the perlin simplex noise
    //


    /**
     * fBM noise.
     *
     * @arg position: The position to seed the noise.
     *
     * @returns: The noise value in the range [-1, 1].
     */
    float fractalBrownianMotionNoise(const float3 &position)
    {
        const float4 position4d = float4(
            position.x,
            position.y,
            position.z,
            0.0f
        );
        float output = 0.0f;
        float frequency = _lacunarity;
        float amplitude = 1.0f;
        float denom = 0.0f;
        float4 translation;
        float4 scale;

        for (int octave=0; octave < _octaves; octave++)
        {
            const float octaveFraction = octave / _octaves;
            scale = (
                (_highFrequencyScale * octaveFraction)
                + (_lowFrequencyScale * (1.0f - octaveFraction))
            );       
            translation = __translation + (
                (_highFrequencyTranslation * octaveFraction)
                + (_lowFrequencyTranslation * (1.0f - octaveFraction))
            );

            output += amplitude * perlinSimplexNoise(
                (position4d * scale + translation) * frequency / _size
            );

            frequency *= _lacunarity;
            denom += amplitude;
            amplitude *= _gain;
        }

        return pow(fabs(output / denom), 1.0f / _gamma);
    }


    /**
     * Turbulence noise.
     *
     * @arg position: The position to seed the noise.
     *
     * @returns: The noise value in the range [0, 1].
     */
    float turbulenceNoise(const float3 &position)
    {
        const float4 position4d = float4(
            position.x,
            position.y,
            position.z,
            0.0f
        );
        float output = 0.0f;
        float frequency = _lacunarity;
        float amplitude = 1.0f;
        float denom = 0.0f;
        float4 translation;
        float4 scale;

        for (int octave=0; octave < _octaves; octave++)
        {
            const float octaveFraction = octave / _octaves;
            scale = (
                (_highFrequencyScale * octaveFraction)
                + (_lowFrequencyScale * (1.0f - octaveFraction))
            );       
            translation = __translation + (
                (_highFrequencyTranslation * octaveFraction)
                + (_lowFrequencyTranslation * (1.0f - octaveFraction))
            );

            output += fabs(
                amplitude * perlinSimplexNoise(
                    (position4d * scale + translation) * frequency / _size
                )
            );

            frequency *= _lacunarity;
            denom += amplitude;
            amplitude *= _gain;
        }

        return pow(output / denom, 1.0f / _gamma);
    }


    /**
     * Create a ray out of the camera. It will be either a standard ray,
     * a latlong ray, or a ray that will result in depth of field.
     *
     * @arg seed: The seed to use in randomization.
     * @arg pixelLocation: The x, and y locations of the pixel.
     * @arg rayOrigin: The location to store the origin of the new ray.
     * @arg rayDirection: The location to store the direction of the new
     *     ray.
     */
    void getCameraRay(
            const float2 &seed,
            const float2 &pixelLocation,
            float3 &rayOrigin,
            float3 &rayDirection)
    {
        const float2 uvCoordinates = pixelsToUV(
            pixelLocation + random(seed),
            float2(_formatWidth, _formatHeight)
        );
        if (_latLong)
        {
            createLatLongCameraRay(
                _cameraWorldMatrix,
                uvCoordinates,
                rayOrigin,
                rayDirection
            );
        }
        else if (_depthOfFieldEnabled)
        {
            createCameraRay(
                _cameraWorldMatrix,
                __inverseCameraProjectionMatrix,
                uvCoordinates,
                __aperture,
                _focalDistance,
                seed,
                rayOrigin,
                rayDirection
            );
        }
        else
        {
            createCameraRay(
                _cameraWorldMatrix,
                __inverseCameraProjectionMatrix,
                uvCoordinates,
                rayOrigin,
                rayDirection
            );
        }
    }


    /**
     * Compute the minimum distance to an object in the scene.
     *
     * @arg rayOrigin: The origin position of the ray.
     * @arg pixelFootprint: A value proportional to the amount of world
     *     space that fills a pixel, like the distance from camera.
     *
     * @returns: The minimum distance to an object in the scene.
     */
    float getMinDistanceToObjectInScene(const float3 &point)
    {
        if (_enableSphericalRamp)
        {
            return distanceToSphere(
                point - _sphericalRampPosition,
                _sphericalRampRadii.y
            );
        }
        else if (_enablePlanarRamp)
        {
            return max(
                distanceToPlane(
                    point - (_planePosition + __planeNormal * min(_planarRamp.x, _planarRamp.w)),
                    -__planeNormal
                ),
                distanceToPlane(
                    point - (_planePosition + __planeNormal * max(_planarRamp.x, _planarRamp.w)),
                    __planeNormal
                )
            );
        }
        else
        {
            const float3 toBoxCenter = point - _boxRampPosition;
            return distanceToRectangularPrism(
                matmul(__inverseBoxRotMatrix, toBoxCenter),
                _boxRampDimensions.x,
                _boxRampDimensions.y,
                _boxRampDimensions.z
            ) - _boxRampFalloff;
        }
    }


    /**
     * Estimate the surface normal at the closest point on the closest
     * object to a point
     *
     * @arg point: The point near which to get the surface normal
     * @arg pixelFootprint: A value proportional to the amount of world
     *     space that fills a pixel, like the distance from camera.
     *
     * @returns: The normalized surface normal.
     */
    float3 estimateSurfaceNormal(const float3 &point)
    {
        return normalize(
            __offset0 * getMinDistanceToObjectInScene(
                point + __offset0 * _hitTolerance
            )
            + __offset1 * getMinDistanceToObjectInScene(
                point + __offset1 * _hitTolerance
            )
            + __offset2 * getMinDistanceToObjectInScene(
                point + __offset2 * _hitTolerance
            )
            + __offset3 * getMinDistanceToObjectInScene(
                point + __offset3 * _hitTolerance
            )
        );
    }


    void computeDepthAndSampleStep(
        const float2 &seed,
        const float3 &rayOrigin,
        const float3 &rayDirection,
        float &initialDepth,
        float &sampleStep
    ) {
        if (_enableDepthRamp)
        {
            sampleStep = __depthSampleStep;
            initialDepth = _depthRamp.x;
        }
        else if (!__enableRaymarching)
        {
            initialDepth = 0.0f;
        }

        if (!__enableRaymarching)
        {
            initialDepth += random(seed.x + seed.y) * sampleStep;
            return;
        }

        float pixelFootprint = _hitTolerance;
        float distance = 0.0f;
        float stepDistance = 0.0f;
        bool hitSurface = false;
        float3 position = rayOrigin;

        while (distance < fabs(_maxDistance))
        {
            const float signedDistance = getMinDistanceToObjectInScene(position);

            if(distance == 0.0f && signedDistance <= 0.0f)
            {
                hitSurface = true;
                initialDepth = 0.0f;
            }

            stepDistance = fabs(signedDistance);
            distance += stepDistance;
            position += stepDistance * rayDirection;

            if (stepDistance < pixelFootprint)
            {
                if (hitSurface)
                {
                    sampleStep = min(
                        sampleStep,
                        (distance - initialDepth) / (float) _samplesPerRay
                    );
                    if (sampleStep == 0.0f)
                    {
                        sampleStep = _maxDistance / (float) _samplesPerRay;
                    }
                    initialDepth += random(seed.x + seed.y) * sampleStep;
                    return;
                }
                hitSurface = true;
                initialDepth = min(initialDepth, distance);

                const float3 normal = estimateSurfaceNormal(position + stepDistance * rayDirection);
                position += 2.0f * pixelFootprint * (rayDirection - normal);
                pixelFootprint = _hitTolerance;
                continue;
            }

            pixelFootprint += _hitTolerance * stepDistance;
        }

        if (hitSurface)
        {
            sampleStep = min(
                sampleStep,
                (distance - initialDepth) / (float) _samplesPerRay
            );
            initialDepth += random(seed.x + seed.y) * sampleStep;
        }
        else
        {
            initialDepth = -1.0f;
        }
    }


    /**
     * Compute a noise value.
     *
     * @arg pos: The x, and y location we are currently processing.
     */
    void process(int2 pos)
    {
        float2 pixelLocation = float2(pos.x, pos.y);
        SampleType(seeds) seedPixel = seeds();
        float2 seed = float2(seedPixel.x, seedPixel.y) + RAND_CONST_0 * pixelLocation;

        float4 resultPixel = float4(0.0f);

        // Generate a ray from the camera
        float3 rayOrigin;
        float3 rayDirection;
        getCameraRay(
            seed,
            pixelLocation,
            rayOrigin,
            rayDirection
        );

        // Set the depth to the start of the depth ramp and add a random offset
        // to eliminate layer lines
        float depth = _maxDistance;
        float sampleStep = _maxDistance / (float) _samplesPerRay;
        computeDepthAndSampleStep(
            seed,
            rayOrigin + rayDirection * _nearPlane,
            rayDirection,
            depth,
            sampleStep
        );

        if (depth < 0.0f)
        {
            dst() = resultPixel;
            return;
        }

        rayOrigin += depth * rayDirection;

        depth += (1.0f + _individualSample) * sampleStep;
        rayOrigin += (1.0f + _individualSample) * rayDirection * sampleStep;
        for (int sample=0; sample <= (int) _secondSample; sample++)
        {
            float ramp = 1.0f;

            float adjustedDepth = depth;
            if (_useCameraDepth)
            {
                adjustedDepth = fabs(matmul(
                    __inverseCameraWorldMatrix,
                    float4(rayOrigin.x, rayOrigin.y, rayOrigin.z, 1.0f)
                )[2]);
            }

            resultPixel[sample * 2 + 1] = adjustedDepth;

            if (_enableSphericalRamp)
            {
                const float minDistanceToSphere = distanceToSphere(
                    rayOrigin - _sphericalRampPosition,
                    _sphericalRampRadii.x
                );
                if (minDistanceToSphere >= __sphericalRampFalloffDistance)
                {
                    resultPixel[sample * 2] = 0.0f;
                    continue;
                }
                else if (minDistanceToSphere > 0.0f)
                {
                    ramp *= (
                        (__sphericalRampFalloffDistance - minDistanceToSphere)
                        / __sphericalRampFalloffDistance
                    );
                }
                if (_sphericalFalloffPower != 0.0f)
                {
                    ramp /= pow(
                        max(
                            1e-6,
                            fabs(
                                _sphericalFalloffOffset
                                + minDistanceToSphere
                                + _sphericalRampRadii.x
                            )
                        ),
                        _sphericalFalloffPower
                    );
                }
            }
            else if (_enablePlanarRamp)
            {
                // Apply the scaling specified by the planar ramp
                const float minDistanceToPlane = distanceToPlane(
                    rayOrigin - _planePosition,
                    __planeNormal
                );
                if (
                    minDistanceToPlane < _planarRamp.x
                    || minDistanceToPlane >= _planarRamp.w
                ) {
                    resultPixel[sample * 2] = 0.0f;
                    continue;
                }
                else if (minDistanceToPlane < _planarRamp.y)
                {
                    ramp *= (
                        (minDistanceToPlane - _planarRamp.y)
                        / (_planarRamp.y - _planarRamp.x)
                        + 1.0f
                    );
                }
                else if (minDistanceToPlane > _planarRamp.z)
                {
                    ramp *= (
                        (_planarRamp.w - minDistanceToPlane)
                        / (_planarRamp.w - _planarRamp.z)
                    );
                }
                if (_planarFalloffPower != 0.0f)
                {
                    ramp /= pow(
                        max(
                            1e-6,
                            fabs(_planarFalloffOffset + minDistanceToPlane)
                        ),
                        _planarFalloffPower
                    );
                }
            }
            else if (_enableBoxRamp)
            {
                const float3 toBoxCenter = rayOrigin - _boxRampPosition;
                const float minDistanceToBox = distanceToRectangularPrism(
                    matmul(__inverseBoxRotMatrix, toBoxCenter),
                    _boxRampDimensions.x,
                    _boxRampDimensions.y,
                    _boxRampDimensions.z
                );
                if (minDistanceToBox >= _boxRampFalloff)
                {
                    resultPixel[sample * 2] = 0.0f;
                    continue;
                }
                else if (minDistanceToBox > 0.0f)
                {
                    ramp *= (_boxRampFalloff - minDistanceToBox) / _boxRampFalloff;
                }
                if (_boxFalloffPower != 0.0f)
                {
                    ramp /= pow(
                        max(
                            1e-6,
                            _boxFalloffOffset + length(toBoxCenter)
                        ),
                        _boxFalloffPower
                    );
                }
            }

            if (_enableDepthRamp)
            {
                // Apply the scaling specified by the depth ramp
                if (adjustedDepth < _depthRamp.y)
                {
                    ramp *= (
                        (adjustedDepth - _depthRamp.y)
                        / (_depthRamp.y - _depthRamp.x)
                        + 1.0f
                    );
                }
                else if (adjustedDepth > _depthRamp.z)
                {
                    ramp *= (_depthRamp.w - adjustedDepth) / (_depthRamp.w - _depthRamp.z);
                }
            }

            // Compute the noise value at this position
            float noiseValue = _density * ramp;
            if (_noiseType == FBM_NOISE)
            {
                noiseValue *= fractalBrownianMotionNoise(rayOrigin);
            }
            else
            {
                noiseValue *= turbulenceNoise(rayOrigin);
            }

            resultPixel[sample * 2] = noiseValue;

            depth += sampleStep;
            rayOrigin += rayDirection * sampleStep;
        }

        dst() = resultPixel;
    }
};
