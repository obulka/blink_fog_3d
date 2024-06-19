// Copyright 2024 by Owen Bulka.
// All rights reserved.
// This file is released under the "MIT License Agreement".
// Please see the LICENSE.md file that should have been included as part
// of this package.


//
// Math operations
//


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



kernel FogKernel : ImageComputationKernel<ePixelWise>
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

        // Image params
        float _formatWidth;
        float _formatHeight;
        float _individualSample;
        float _density;
        int _samplesPerRay;
        float4 _depthRamp;
        float4 _yRamp;
        bool _enableYRamp;
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

    local:
        // These local variables are not exposed to the user.
        float4x4 __inverseCameraProjectionMatrix;
        float __aperture;

        int __simplex[64][4];
        int __perm[512];
        int __grad4[32][4];

    /**
     * Give the parameters labels and default values.
     */
    void define()
    {
        // Camera params
        defineParam(_focalLength, "Focal Length", 50.0f);
        defineParam(_horizontalAperture, "Horizontal Aperture", 24.576f);
        defineParam(_nearPlane, "Near Plane", 0.1f);
        defineParam(_farPlane, "Far Plane", 10000.0f);
        defineParam(
            _cameraWorldMatrix,
            "Camera World Matrix",
            float4x4(
                1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1
            )
        );
        defineParam(_focalDistance, "Focal Distance", 4.0f);
        defineParam(_fStop, "fstop", 16.0f);
        defineParam(_depthOfFieldEnabled, "Enable Depth Of Field", true);
        defineParam(_latLong, "Output LatLong", false);

        // Image params
        defineParam(_formatHeight, "Screen Height", 2160.0f);
        defineParam(_formatWidth, "Screen Width", 3840.0f);
        defineParam(_individualSample, "Individual Sample", 0.0f);
        defineParam(_density, "Density", 1.0f);
        defineParam(_samplesPerRay, "Samples Per Ray", 5);
        defineParam(_depthRamp, "Sample Depth Ramp", float4(5.0f, 10.0f, 15.0f, 20.0f));
        defineParam(_yRamp, "Sample Y Ramp", float4(5.0f, 10.0f, 15.0f, 20.0f));
        defineParam(_enableYRamp, "Enable Y Ramp", false);
        defineParam(_secondSample, "Do Second Sample", false);

        // Noise Parameters
        defineParam(_size, "Size", 20.0f);
        defineParam(_noiseType, "Noise Type", 0);
        defineParam(_translation, "Translation", float3(0.0f, 0.0f, 0.0f));
        defineParam(_octaves, "Octaves", 8.0f);
        defineParam(_lacunarity, "Lacunarity", 3.0f);
        defineParam(_gain, "Gain", 0.5f);
        defineParam(_gamma, "Gamma", 0.5f);
        defineParam(_lowFrequencyScale, "Low Frequency Scale", float4(0.0f, 0.0f, 0.0f, 0.0f));
        defineParam(_highFrequencyScale, "High Frequency Scale", float4(0.0f, 0.0f, 0.0f, 0.0f));
        defineParam(_lowFrequencyTranslation, "Low Frequency Translation", float4(0.0f, 0.0f, 0.0f, 0.0f));
        defineParam(_highFrequencyTranslation, "High Frequency Translation", float4(0.0f, 0.0f, 0.0f, 0.0f));
    }


    /**
     * Initialize the local variables.
     */
    void init()
    {
        const float aspect = aspectRatio(_formatHeight, _formatWidth);
        float4x4 cameraProjectionMatrix = projectionMatrix(
            _focalLength,
            _horizontalAperture,
            aspect,
            _nearPlane,
            _farPlane
        );
        __inverseCameraProjectionMatrix = cameraProjectionMatrix.invert();

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
     * @arg simplex: The simplex LUT.
     * @arg perm: The perm LUT.
     * @arg grad4: The grad4 LUT.
     *
     * @returns: Noise value in the range [-1, 1], value of 0 on all integer
     *     coordinates.
     */
    inline float perlinSimplexNoise(
            const float4 &seed
    )
    {
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
    float fractalBrownianMotionNoise(const float4 &position)
    {
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
            translation = (
                (_highFrequencyTranslation * octaveFraction)
                + (_lowFrequencyTranslation * (1.0f - octaveFraction))
            );

            output += amplitude * perlinSimplexNoise(
                (position * scale + translation) * frequency / _size
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
    float turbulenceNoise(const float4 &position)
    {
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
            translation = (
                (_highFrequencyTranslation * octaveFraction)
                + (_lowFrequencyTranslation * (1.0f - octaveFraction))
            );

            output += fabs(
                amplitude * perlinSimplexNoise(
                    (position * scale + translation) * frequency / _size
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

        const float sampleStep = (_depthRamp.w - _depthRamp.x) / (float) _samplesPerRay;

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
        float depth = _depthRamp.x - random(seed.x + seed.y + _individualSample) * sampleStep;
        rayOrigin += depth * rayDirection + _translation;

        depth += (1.0f + _individualSample) * sampleStep;
        rayOrigin += (1.0f + _individualSample) * rayDirection * sampleStep;
        for (int sample=0; sample <= (int) _secondSample; sample++)
        {
            resultPixel[sample * 2 + 1] = depth;

            // Get the noise value based on which type of noise we are using
            const float4 samplePosition4d = float4(
                rayOrigin.x,
                rayOrigin.y,
                rayOrigin.z,
                0.0f
            );

            // Apply the scaling specified by the y ramp
            float yRamp = 1.0f;
            if (_enableYRamp)
            {
                const float yPosition = samplePosition4d.y;
                if (yPosition <= _yRamp.x || yPosition >= _yRamp.w)
                {
                    resultPixel[sample * 2] = 0.0f;
                    continue;
                }
                else if (yPosition < _yRamp.y)
                {
                    yRamp = (yPosition - _yRamp.y) / (_yRamp.y - _yRamp.x) + 1.0f;
                }
                else if (yPosition > _yRamp.z)
                {
                    yRamp = (_yRamp.w - yPosition) / (_yRamp.w - _yRamp.z);
                }
            }

            // Apply the scaling specified by the depth ramp
            float depthRamp = 1.0f;
            if (depth < _depthRamp.y)
            {
                depthRamp = (depth - _depthRamp.y) / (_depthRamp.y - _depthRamp.x) + 1.0f;
            }
            else if (depth > _depthRamp.z)
            {
                depthRamp = (_depthRamp.w - depth) / (_depthRamp.w - _depthRamp.z);
            }

            // Compute the noise value at this position
            float noiseValue = _density * depthRamp * yRamp;
            if (_noiseType == FBM_NOISE)
            {
                noiseValue *= fractalBrownianMotionNoise(samplePosition4d);
            }
            else
            {
                noiseValue *= turbulenceNoise(samplePosition4d);
            }

            resultPixel[sample * 2] = noiseValue;

            depth += sampleStep;
            rayOrigin += rayDirection * sampleStep;
        }

        dst() = resultPixel;
    }
};
