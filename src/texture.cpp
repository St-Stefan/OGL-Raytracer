#include "texture.h"
#include <framework/image.h>
#include <glm/geometric.hpp>
#include <shading.h>


// Source: Chapter 11.1, Fundamentals of Computer Graphics
glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    int i = std::roundf(texCoord.x * image.width - 0.5f);
    int j = std::roundf((1 - texCoord.y) * image.height - 0.5f);
    int index = j * image.width + i;
    return acquireTexelClampMode(i, j, image);
}

glm::vec3 acquireTexelClampMode(int i, int j, const Image& image)
{
    if (j * image.width + i >= image.pixels.size()) {
        return image.pixels[image.pixels.size() - 1];
    }
    if (j * image.width + i < 0) {
        return image.pixels[0];
    }
    return image.pixels[j * image.width + i];
}

std::pair<std::pair<float,float>, int> cubeMapLookUp(float x, float y, float z)
{
    float absX = fabs(x);
    float absY = fabs(y);
    float absZ = fabs(z);
    int index = 0; 
    float u = 0.0f; 
    float v = 0.0f;

    int isXPositive = x > 0 ? 1 : 0;
    int isYPositive = y > 0 ? 1 : 0;
    int isZPositive = z > 0 ? 1 : 0;

    float maxAxis = 0.0f;
    float uc = 0.0f;
    float vc = 0.0f;

    // POSITIVE X
    if (isXPositive && absX >= absY && absX >= absZ) {
        // u (0 to 1) goes from +z to -z
        // v (0 to 1) goes from -y to +y
        maxAxis = absX;
        uc = -z;
        vc = y;
        index = 0;
    }
    // NEGATIVE X
    if (!isXPositive && absX >= absY && absX >= absZ) {
        // u (0 to 1) goes from -z to +z
        // v (0 to 1) goes from -y to +y
        maxAxis = absX;
        uc = z;
        vc = y;
        index = 1;
    }
    // POSITIVE Y
    if (isYPositive && absY >= absX && absY >= absZ) {
        // u (0 to 1) goes from -x to +x
        // v (0 to 1) goes from +z to -z
        maxAxis = absY;
        uc = x;
        vc = -z;
        index = 2;
    }
    // NEGATIVE Y
    if (!isYPositive && absY >= absX && absY >= absZ) {
        // u (0 to 1) goes from -x to +x
        // v (0 to 1) goes from -z to +z
        maxAxis = absY;
        uc = x;
        vc = z;
        index = 3;
    }
    // POSITIVE Z
    if (isZPositive && absZ >= absX && absZ >= absY) {
        // u (0 to 1) goes from -x to +x
        // v (0 to 1) goes from -y to +y
        maxAxis = absZ;
        uc = x;
        vc = y;
        index = 4;
    }
    // NEGATIVE Z
    if (!isZPositive && absZ >= absX && absZ >= absY) {
        // u (0 to 1) goes from +x to -x
        // v (0 to 1) goes from -y to +y
        maxAxis = absZ;
        uc = -x;
        vc = y;
        index = 5;
    }

    // Convert range from -1 to 1 to 0 to 1
    u = 0.5f * (uc / maxAxis + 1.0f);
    v = 0.5f * (vc / maxAxis + 1.0f);

    return {{ u, v }, index};
}

glm::vec3 environmentMapping(const std::vector<Image>& images,
        const Ray& ray, const Features& features)
{  
    std::pair<std::pair<float, float>, int> result = cubeMapLookUp(ray.direction.x, ray.direction.y, ray.direction.z);
    return acquireTexel(images[result.second], glm::vec2{result.first.first, result.first.second} , features);
}